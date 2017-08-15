#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_first_of("}");
    if (found_null != string::npos) {
        return "";
    } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y) {

    double closestLen = 100000; //large number
    int closestWaypoint = 0;

    for (int i = 0; i < maps_x.size(); i++) {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x, y, map_x, map_y);
        if (dist < closestLen) {
            closestLen = dist;
            closestWaypoint = i;
        }

    }

    return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {

    int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];

    double heading = atan2((map_y - y), (map_x - x));

    double angle = abs(theta - heading);

    if (angle > pi() / 4) {
        closestWaypoint++;
    }

    return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {
    int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

    int prev_wp;
    prev_wp = next_wp - 1;
    if (next_wp == 0) {
        prev_wp = maps_x.size() - 1;
    }

    double n_x = maps_x[next_wp] - maps_x[prev_wp];
    double n_y = maps_y[next_wp] - maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
    double proj_x = proj_norm * n_x;
    double proj_y = proj_norm * n_y;

    double frenet_d = distance(x_x, x_y, proj_x, proj_y);

    //see if d value is positive or negative by comparing it to a center point

    double center_x = 1000 - maps_x[prev_wp];
    double center_y = 2000 - maps_y[prev_wp];
    double centerToPos = distance(center_x, center_y, x_x, x_y);
    double centerToRef = distance(center_x, center_y, proj_x, proj_y);

    if (centerToPos <= centerToRef) {
        frenet_d *= -1;
    }

    // calculate s value
    double frenet_s = 0;
    for (int i = 0; i < prev_wp; i++) {
        frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
    }

    frenet_s += distance(0, 0, proj_x, proj_y);

    return {frenet_s, frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y) {
    int prev_wp = -1;

    while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
        prev_wp++;
    }

    int wp2 = (prev_wp + 1) % maps_x.size();

    double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s - maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
    double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

    double perp_heading = heading - pi() / 2;

    double x = seg_x + d * cos(perp_heading);
    double y = seg_y + d * sin(perp_heading);

    return {x, y};

}

void checkTrafficAhead(int ego_car_lane, double ego_car_s, double time_horizon,
                       const vector<vector<double>> &sensor_fusion, bool *too_close, double *car_ahead_speed) {
    for (int i = 0; i < sensor_fusion.size(); ++i) {
        float d = sensor_fusion[i][6];
        double check_car_s = sensor_fusion[i][5];
        if (d < (2 + ego_car_lane * 4 + 2) && d > (2 + ego_car_lane * 4 - 2)) {
            double vx = sensor_fusion[i][3];
            double vy = sensor_fusion[i][4];
            double check_car_speed = sqrt(pow(vx, 2) + pow(vy, 2));
            check_car_s += time_horizon * check_car_speed;
            if (check_car_s > ego_car_s && check_car_s - ego_car_s < 20) {
                *too_close = true;
                *car_ahead_speed = check_car_speed;
                break;
            }
        }
    }
}

void findNewLane(int ego_car_lane, double ego_car_s, double time_horizon, const vector<vector<double>> &sensor_fusion,
                 int *new_lane) {
    for (int l = 0; l <= 2; ++l) {
        if ((l != ego_car_lane) && abs(l - ego_car_lane) < 2) {
            bool safe_to_change_lanes = true;
            for (int i = 0; i < sensor_fusion.size(); ++i) {
                float d = sensor_fusion[i][6];
                if (d < (2 + l * 4 + 2) && d > (2 + l * 4 - 2)) {
                    double vx = sensor_fusion[i][3];
                    double vy = sensor_fusion[i][4];
                    double speed = sqrt(pow(vx, 2) + pow(vy, 2));
                    double check_car_s = sensor_fusion[i][5];
                    check_car_s += time_horizon * speed;
                    if (abs(check_car_s - ego_car_s) < 20) {
                        safe_to_change_lanes = false;
                        break;
                    }
                }
            }
            if (safe_to_change_lanes) {
                *new_lane = l;
                break;
            }
        }
    }
}

void addPreviousPathPoints(const vector<double> &previous_path_x, const vector<double> &previous_path_y,
                           vector<double> *next_x_vals, vector<double> *next_y_vals) {
    for (int i = 0; i < previous_path_x.size(); ++i) {
        next_x_vals->push_back(previous_path_x[i]);
        next_y_vals->push_back(previous_path_y[i]);
    }
}

void getAnchorPoints(double ego_car_x, double ego_car_y, double ego_car_yaw, double ego_car_s, double ego_car_lane,
                     const vector<double> &previous_path_x, const vector<double> &previous_path_y,
                     const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y,
                     vector<double> *ptsx, vector<double> *ptsy) {
    int prev_size = previous_path_x.size();
    if (prev_size < 1) {
        ptsx->push_back(ego_car_x - cos(ego_car_yaw));
        ptsy->push_back(ego_car_y - sin(ego_car_yaw));
        ptsx->push_back(ego_car_x);
        ptsy->push_back(ego_car_y);
    } else {
        ptsx->push_back(previous_path_x[prev_size - 2]);
        ptsy->push_back(previous_path_y[prev_size - 2]);
        ptsx->push_back(previous_path_x[prev_size - 1]);
        ptsy->push_back(previous_path_y[prev_size - 1]);
    }

    vector<double> next_wp0 = getXY(ego_car_s + 30, (2 + 4 * ego_car_lane), maps_s, maps_x, maps_y);
    vector<double> next_wp1 = getXY(ego_car_s + 60, (2 + 4 * ego_car_lane), maps_s, maps_x, maps_y);
    vector<double> next_wp2 = getXY(ego_car_s + 90, (2 + 4 * ego_car_lane), maps_s, maps_x, maps_y);

    ptsx->push_back(next_wp0[0]);
    ptsx->push_back(next_wp1[0]);
    ptsx->push_back(next_wp2[0]);

    ptsy->push_back(next_wp0[1]);
    ptsy->push_back(next_wp1[1]);
    ptsy->push_back(next_wp2[1]);
}

vector<double> getCarPose(const vector<double> &prev_path_x, const vector<double> &prev_path_y) {
    int prev_size = prev_path_x.size();
    double x = prev_path_x[prev_size - 1];
    double y = prev_path_y[prev_size - 1];
    double x_prev = prev_path_x[prev_size - 2];
    double y_prev = prev_path_y[prev_size - 2];
    double yaw = atan2(y - y_prev, x - x_prev);
    return {x, y, yaw};
}

void transformAnchorPoints(vector<double> *ptsx, vector<double> *ptsy, double ref_x, double ref_y, double ref_yaw) {
    // Transform spline anchor points to vehicle reference frame
    for (int i = 0; i < ptsx->size(); ++i) {
        double shift_x = (*ptsx)[i] - ref_x;
        double shift_y = (*ptsy)[i] - ref_y;
        (*ptsx)[i] = shift_x * cos(ref_yaw) + shift_y * sin(ref_yaw);
        (*ptsy)[i] = -1 * shift_x * sin(ref_yaw) + shift_y * cos(ref_yaw);
    }
}

void interpolateSpline(double ego_car_target_vel, double ref_x, double ref_y, double ref_yaw, const tk::spline &s,
                       const vector<double> &previous_path_x, vector<double> *next_x_vals,
                       vector<double> *next_y_vals) {
    double target_x = 30.0;
    double target_y = s(target_x);
    double target_dist = sqrt(pow(target_x, 2) + pow(target_y, 2));
    double x_add_on = 0.0;
    double N = (target_dist / (0.02 * ego_car_target_vel * 0.447));
    double world_x_point, world_y_point;
    for (int i = 1; i <= 50 - previous_path_x.size(); i++) {
        double x_point = x_add_on + (target_x / N);
        double y_point = s(x_point);
        x_add_on = x_point;

        // Transform spline point back to world coordinate frame.

        world_x_point = x_point * cos(ref_yaw) - y_point * sin(ref_yaw) + ref_x;
        world_y_point = x_point * sin(ref_yaw) + y_point * cos(ref_yaw) + ref_y;

        // Add spline point to path.
        next_x_vals->push_back(world_x_point);
        next_y_vals->push_back(world_y_point);
    }
}

int main() {
    uWS::Hub h;

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;

    // Waypoint map to read from
    string map_file_ = "../data/highway_map.csv";
    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;

    ifstream in_map_(map_file_.c_str(), ifstream::in);

    string line;
    while (getline(in_map_, line)) {
        istringstream iss(line);
        double x;
        double y;
        float s;
        float d_x;
        float d_y;
        iss >> x;
        iss >> y;
        iss >> s;
        iss >> d_x;
        iss >> d_y;
        map_waypoints_x.push_back(x);
        map_waypoints_y.push_back(y);
        map_waypoints_s.push_back(s);
        map_waypoints_dx.push_back(d_x);
        map_waypoints_dy.push_back(d_y);
    }

    // Lane occupied by vehicle.
    // Lane 0 is the most disant from highway dividing line (rightmost lane).
    // Lane 1 is the middle lane.
    // Lane 2 is the lane closest to the dividing line (leftmost lane).
    int lane = 1;

    // Vehicle target velocity
    double target_vel = 0; // mph

    // Vehicle state
    string state = "KL";


    h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy,
                        &lane, &target_vel, &state](
            uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
            uWS::OpCode opCode) {



        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        //auto sdata = string(data).substr(0, length);
        //cout << sdata << endl;
        if (length && length > 2 && data[0] == '4' && data[1] == '2') {

            auto s = hasData(data);

            if (s != "") {
                auto j = json::parse(s);

                string event = j[0].get<string>();

                if (event == "telemetry") {
                    // j[1] is the data JSON object

                    // Main car's localization Data
                    double car_x = j[1]["x"];
                    double car_y = j[1]["y"];
                    double car_s = j[1]["s"];
                    double car_d = j[1]["d"];
                    double car_yaw = j[1]["yaw"];
                    double car_speed = j[1]["speed"];

                    // Previous path data given to the Planner
                    auto previous_path_x = j[1]["previous_path_x"];
                    auto previous_path_y = j[1]["previous_path_y"];
                    // Previous path's end s and d values
                    double end_path_s = j[1]["end_path_s"];
                    double end_path_d = j[1]["end_path_d"];

                    // Sensor Fusion Data, a list of all other cars on the same side of the road.
                    auto sensor_fusion = j[1]["sensor_fusion"];

                    // Set car_s to be at end of previous path
                    int prev_size = previous_path_x.size();
                    if (prev_size > 0) {
                        car_s = end_path_s;
                    }

                    // Check for vehicles ahead of us that are too close (within 20 meters).
                    bool too_close = false;
                    double car_ahead_speed;
                    double time_horizon = 0.02 * prev_size;
                    checkTrafficAhead(lane, car_s, time_horizon, sensor_fusion, &too_close, &car_ahead_speed);

                    // If a vehicle is ahead of us and too close, then decelerate until we match its speed.
                    if (too_close) {
                        if (target_vel * 0.447 > car_ahead_speed) {
                            target_vel -= 0.30;
                        }
                    }

                    // If no vehicle ahead and in KL state, then accelerate to the target speed of 49.5 mph.
                    if (!too_close) {
                        if (target_vel < 49.5) {
                            target_vel += 0.30;
                        }
                    }

                    // If current speed is below target speed of 49.5 mph then transition from KL state to CL state.
                    if (too_close && target_vel < 49.5 && state.compare("KL") == 0) {
                        state = "CL";
                        findNewLane(lane, car_s, time_horizon, sensor_fusion, &lane);
                    }

                    // If lane change complete, then return to keep-lane state
                    if (state.compare("CL") == 0 && abs(2 + 4 * lane - car_d) < 0.5) {
                        state = "KL";
                    }

                    // Begin new path using points from previous path
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;
                    addPreviousPathPoints(previous_path_x, previous_path_y, &next_x_vals, &next_y_vals);

                    // Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m.
                    // Later we'll interpolate between these points
                    vector<double> ptsx;
                    vector<double> ptsy;

                    // Get Spline Anchor Points
                    getAnchorPoints(car_x, car_y, car_yaw, car_s, lane, previous_path_x, previous_path_y,
                                    map_waypoints_s, map_waypoints_x, map_waypoints_y, &ptsx, &ptsy);


                    // Transform anchor points to vehicle coordinate frame
                    double ref_x, ref_y, ref_yaw;
                    if (prev_size < 2) {
                        ref_x = car_x;
                        ref_y = car_y;
                        ref_yaw = car_yaw;
                    } else {
                        vector<double> pose = getCarPose(previous_path_x, previous_path_y);
                        ref_x = pose[0];
                        ref_y = pose[1];
                        ref_yaw = pose[2];
                    }
                    transformAnchorPoints(&ptsx, &ptsy, ref_x, ref_y, ref_yaw);

                    // Assemble Spline
                    tk::spline s;
                    s.set_points(ptsx, ptsy);

                    // Calculate Spline points to travel at desired velocity given 0.02 sec between successive waypoints
                    interpolateSpline(target_vel, ref_x, ref_y, ref_yaw, s, previous_path_x, &next_x_vals,
                                      &next_y_vals);

                    json msgJson;

                    // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;

                    auto msg = "42[\"control\"," + msgJson.dump() + "]";

                    //this_thread::sleep_for(chrono::milliseconds(1000));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

                }
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });

    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                       size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
















































































