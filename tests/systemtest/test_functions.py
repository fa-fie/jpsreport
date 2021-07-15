﻿import os
import logging
import math
import numpy as np
from utils import FAILURE

###### GENERAL FUNCTIONS #############################

def get_num_time_intervals(num_frames, diff_frames):
    if (num_frames / diff_frames).is_integer():
        return int(num_frames / diff_frames - 1)
    else:
        return int(num_frames / diff_frames)

def get_velocity_range(distance, real_v, fps, tolerance):
    # explanation for accepted range:
    # calculation of number of frames which is needed to pass the length of a measurement area (real_frames)
    # if this value is not an integer, the number of frames in which the pedestrian passes the measurement
    # area could differ (depending on starting position)
    # counted frames should either be the next lower integer or the next higher integer value
    # depending on these frame values, the velocity range is calculated

    real_frames = distance / real_v * fps
    if real_frames.is_integer():
        return [real_v, real_v]
    elif math.floor(real_frames) == 0:
        # combination of these factors: distance too small, real_v too large, fps too small
        logging.critical("jpsreport might not detect velocity correctly; minimum of detected frames is 0")
        max_seconds = math.ceil(real_frames) / fps
        return [distance/max_seconds, distance/max_seconds]
    else:
        min_seconds = math.floor(real_frames) / fps
        max_seconds = math.ceil(real_frames) / fps
        return [distance/max_seconds - tolerance, distance/min_seconds + tolerance]

def get_density_range(distance, real_v, fps, tolerance, num_peds, delta_t):
    real_frames = distance / real_v * fps
    if real_frames.is_integer():
        density = ((num_peds * real_frames) / fps) / (delta_t * distance)
        return [density, density]
    else:
        min_seconds = math.floor(real_frames) / fps
        max_seconds = math.ceil(real_frames) / fps
        return [(num_peds * min_seconds) / (delta_t * distance) - tolerance, 
                (num_peds * max_seconds) / (delta_t * distance) + tolerance]

    if math.floor(real_frames) == 0:
        # combination of these factors: distance too small, real_v too large, fps too small
        logging.critical("jpsreport might not detect density correctly; minimum of detected frames is 0")

def get_num_pass_lines(line_pos, 
                       start_pos_x, dist, num_columns, peds_y, velocity,
                       num_frames, delta_t_frames, delta_t_seconds):
    number_time_intervals = get_num_time_intervals(num_frames, delta_t_frames)
    fps = delta_t_frames/delta_t_seconds
    dist_per_delta_t = velocity * delta_t_frames / fps

    num_peds_per_line = []
    for x_line in line_pos:
        column_position_t0 = [start_pos_x - num_column * dist for num_column in range(num_columns)]
        column_position_t1 = [start_pos_x - num_column * dist + dist_per_delta_t for num_column in range(num_columns)]

        num_columns_delta_t = []
        for i in range(number_time_intervals):
            current_num_columns = 0
            for column_idx in range(num_columns):
                x_pos_t0 = column_position_t0[column_idx]
                x_pos_t1 = column_position_t1[column_idx]

                if x_pos_t0 <= x_line and x_pos_t1 >= x_line:
                    current_num_columns += 1

            column_position_t0 = column_position_t1
            column_position_t1 = [pos + dist_per_delta_t for pos in column_position_t0]
            num_columns_delta_t.append(current_num_columns * peds_y)
        num_peds_per_line.append(num_columns_delta_t)

    return num_peds_per_line

###### FUNCTIONS FOR METHOD E ########################

def get_num_in_area(start_pos_x, dist, num_columns, peds_y, velocity, fps, num_frames, x0, x1):
    dist_per_frame = velocity / fps
    column_position = [start_pos_x - num_column * dist for num_column in range(num_columns)]

    num_peds_per_frame = []
    for i in range(num_frames):
        num_columns_frame = 0
        for column_idx in range(num_columns):
            x_pos = column_position[column_idx]
            scenarios = [
                math.isclose(x_pos, x0, abs_tol=0.00001) or math.isclose(x_pos, x1, abs_tol=0.00001),
                x_pos >= x0 and x_pos <= x1,
                x_pos > x1 and x_pos < (x1 + dist_per_frame) and not math.isclose(x_pos - dist_per_frame, x1, abs_tol=0.00001)
                ]
            # last condition is added because of variant used to detect tIn/tOut
            # this has to be modified if another variant is used
            if True in scenarios:
                num_columns_frame += 1
        tmp_pos = column_position
        column_position = [pos + dist_per_frame for pos in tmp_pos]
        num_peds_per_frame.append(num_columns_frame * peds_y)

    return num_peds_per_frame

def check_velocity_values_E(out_v, general_reference_rho, reference_specific_flow, delta_t_frames, 
                            num_frame_intervals, abs_tol):
    rho_time_interval_idx = []
    if num_frame_intervals > 1:
        frame_nr = 0
        reference_rho = []
        for val in general_reference_rho:
            reference_rho.append(val)
            if frame_nr % delta_t_frames == 0 and frame_nr != 0:
                reference_rho.append(val)
                rho_time_interval_idx.append(frame_nr)
            frame_nr += 1
        rho_time_interval_idx.pop()
    else:
        reference_rho = general_reference_rho

    success = True
    idx_rho = 0
    time_interval_idx = 0
    for v_value in out_v:
        if not (reference_rho[idx_rho] == 0 or reference_specific_flow[time_interval_idx] == 0):
            ref_v = reference_specific_flow[time_interval_idx] / reference_rho[idx_rho]
        elif reference_rho[idx_rho] == 0 and reference_specific_flow[time_interval_idx] == 0:
            ref_v = math.nan
        elif reference_rho[idx_rho] == 0:
            ref_v = math.inf
        success = success and (math.isclose(ref_v, v_value, abs_tol=abs_tol) or (np.isnan(ref_v) and np.isnan(v_value)))
        # isclose returns True if ref_v and v_value are inf, but not if they are nan

        if idx_rho in rho_time_interval_idx:
            time_interval_idx += 1
        idx_rho += 1
    return success

def runtest_method_E(trajfile,
                     delta_t_frames, delta_t_seconds, num_frames,
                     delta_x, delta_y, line_ids, 
                     real_velocity, 
                     fps,
                     number_pass_line, number_in_area):
    success = True
    abs_tolerance = 0.0001
    # values are rounded differently in output files -> +- 0.0001 as tolerance
    general_output_path = os.path.join('./Output', 'Fundamental_Diagram', 'Method_E')

    logging.info("===== Method E =========================")

    # check density values separate from other values, as those are contained in only one file

    out_rho_fname = os.path.join(general_output_path, f'rho_{trajfile}_id_1.dat')
    if not os.path.exists(out_rho_fname):
        # check if file was output correctly
        logging.critical("jpsreport did not output results correctly.")
        exit(FAILURE)

    out_rho = np.genfromtxt(out_rho_fname)
    reference_rho_dx = np.array([num_peds / delta_x for num_peds in number_in_area])
    reference_rho_area = np.array([num_peds / (delta_x * delta_y) for num_peds in number_in_area])
    if not out_rho.shape == (num_frames, 3):
        success = False
        logging.critical('wrong number of frames for density values')
    else:
        logging.info('correct number of frames for density values')

        if not np.allclose(reference_rho_dx, out_rho[:, 1], atol=abs_tolerance):
            success = False
            logging.critical('wrong density values (delta x)')
        else:
            logging.info('correct density values (delta x)')
        if not np.allclose(reference_rho_area, out_rho[:, 2], atol=abs_tolerance):
            success = False
            logging.critical('wrong density values (area)')
        else:
            logging.info('correct density values (area)')

    # check flow and velocity values per line
    reference_flow_lines = [[num_peds / delta_t_seconds for num_peds in line] for line in number_pass_line]
    reference_specific_flow_lines = [[num_peds / (delta_t_seconds * delta_y) for num_peds in line] for line in number_pass_line]

    for line_id in line_ids:
        logging.info(f'---- line id {line_id} ----')

        out_flow_fname = os.path.join(general_output_path, 
                                      f'flow_{trajfile}_id_1_line_{line_id}.dat')
        out_v_fname = os.path.join(general_output_path, 
                                   f'v_{trajfile}_id_1_line_{line_id}.dat')

        if not os.path.exists(out_flow_fname) or not os.path.exists(out_v_fname):
            # check if files were output correctly
            logging.critical("jpsreport did not output results correctly.")
            exit(FAILURE)

        out_flow = np.genfromtxt(out_flow_fname)
        out_v = np.genfromtxt(out_v_fname)

        number_time_intervals = get_num_time_intervals(num_frames, delta_t_frames)

        # check flow values
        if not (out_flow.shape == (number_time_intervals, 2) and number_time_intervals > 1) and \
            not (out_flow.shape == (2,) and number_time_intervals == 1):
            success = False
            logging.critical('wrong number of time intervals for flow values')
        else:
            logging.info('correct number of time intervals for flow values')

            reference_flow = np.array(reference_flow_lines[line_ids.index(line_id)])
            reference_specific_flow = np.array(reference_specific_flow_lines[line_ids.index(line_id)])
            if number_time_intervals == 1:
                condition_flow = not math.isclose(reference_flow[0], out_flow[0], abs_tol=abs_tolerance)
                condition_specific_flow = not math.isclose(reference_specific_flow[0], out_flow[1], abs_tol=abs_tolerance)
            else:
                condition_flow = not np.allclose(reference_flow, out_flow[:, 0], atol=abs_tolerance)
                condition_specific_flow = not np.allclose(reference_specific_flow, out_flow[:, 1], atol=abs_tolerance)
            
            if condition_flow:
                success = False
                logging.critical('wrong flow values')
            else:
                logging.info('correct flow values')
            if condition_specific_flow:
                success = False
                logging.critical('wrong specific flow values')
            else:
                logging.info('correct specific flow values')

            # check velocity values
            number_velocity_values = number_time_intervals * (delta_t_frames + 1)

            # for multiple delta t -> must check velocity differently!

            if not out_v.shape == (number_velocity_values, 2):
                success = False
                logging.critical('wrong number of frames for velocity values')
            else:
                logging.info('correct number of frames for velocity values')

                check_v = check_velocity_values_E(out_v[:, 1], reference_rho_area, reference_specific_flow, 
                                                  delta_t_frames, number_time_intervals, abs_tolerance)
                if not check_v:
                    success = False
                    logging.critical('wrong velocity values')
                else:
                    logging.info('correct velocity values')

    return success

###### FUNCTIONS FOR METHOD F ########################

def check_velocity_values_F(out_v_values, accepted_range, abs_tolerance):
    success = True
    for v_value in out_v_values:
        if accepted_range[1] == accepted_range[0]:
            condition = not math.isclose(v_value, accepted_range[0], abs_tol=abs_tolerance)
        else:
            condition = accepted_range[0] > v_value or v_value > accepted_range[1]

        if condition:
            success = False
            logging.critical(f'velocity value ({v_value}) is not in accepted range around {real_velocity}m/s')
    return success

def runtest_method_F(trajfile,
                     delta_t_frames, delta_t_seconds, num_frames,
                     delta_x, delta_y, line_ids, 
                     real_velocity, 
                     fps,
                     ped_ids, number_pass_line):
    success = True
    abs_tolerance = 0.0001 
    # values are rounded differently in output files -> +- 0.0001 as tolerance

    accepted_range = get_velocity_range(delta_x, real_velocity, fps, abs_tolerance)
    # range in which velocity values are considered as correct

    general_output_path = os.path.join('./Output', 'Fundamental_Diagram', 'Method_F')
    
    logging.info("===== Method F =========================")

    # check velocity file
    out_v_fname = os.path.join(general_output_path, f'v_{trajfile}_id_1.dat')
    if not os.path.exists(out_v_fname):
        logging.critical("jpsreport did not output results correctly.")
        exit(FAILURE)
    out_v = np.loadtxt(out_v_fname)

    # check velocity values for validity
    if out_v.shape != (2,):
        vel_check = check_velocity_values_F(out_v[:, 1], accepted_range, abs_tolerance)
        out_v_ids = list(out_v[:, 0])
    else:
        vel_check = check_velocity_values_F([out_v[1]], accepted_range, abs_tolerance)
        out_v_ids = [out_v[0]]

    if vel_check:
        logging.info('correct velocity values')

    # check velocity pedestrian ids for validity
    if sorted(out_v_ids) != ped_ids:
        success = False
        logging.critical('wrong pedestrian ids for velocity were output')
    else:
        logging.info('correct pedestrian ids for velocity')

    # check density values per line
    for line_id in line_ids:
        out_rho_flow_fname = os.path.join(general_output_path, 
                                            f'rho_flow_{trajfile}_id_1_line_{line_id}.dat')

        if not os.path.exists(out_rho_flow_fname):
            logging.critical("jpsreport did not output results correctly.")
            exit(FAILURE)

        out_rho_flow = np.loadtxt(out_rho_flow_fname)

        number_time_intervals = get_num_time_intervals(num_frames, delta_t_frames)

        if not (out_rho_flow.shape == (number_time_intervals, 4) and number_time_intervals > 1) and \
            not (out_rho_flow.shape == (4,) and number_time_intervals == 1):
            success = False
            logging.critical(f'line id <{line_id}> wrong number of time intervals')
        else:
            logging.info(f'line id <{line_id}> correct number of time intervals')

            if number_time_intervals == 1:
                out_num_peds = [out_rho_flow[0]]
                out_rho = [out_rho_flow[1]]
                out_flow = [out_rho_flow[2]]
                out_specific_flow = [out_rho_flow[3]]
            else:
                out_num_peds = list(out_rho_flow[:, 0])
                out_rho = out_rho_flow[:, 1]
                out_flow = out_rho_flow[:, 2]
                out_specific_flow = out_rho_flow[:, 3]

            # check number of pedestrians
            if out_num_peds != number_pass_line[line_ids.index(line_id)]:
                success = False
                logging.critical(f'line id <{line_id}> different number of pedestrians counted'
                                 f'({out_num_peds}, not {number_pass_line[line_ids.index(line_id)]})')
            else:
                logging.info(f'line id <{line_id}> correct number of pedestrians counted')
                
            # check flow value
            ref_flow = np.array([num_peds_delta_t / delta_t_seconds for num_peds_delta_t in 
                                 number_pass_line[line_ids.index(line_id)]])
            if not np.allclose(ref_flow, out_flow, atol=abs_tolerance):
                success = False
                logging.critical(f'line id <{line_id}> wrong flow value')
            else:
                logging.info(f'line id <{line_id}> correct flow value')
            
            # check specific flow value
            ref_specific_flow = np.array([num_peds_delta_t / (delta_t_seconds * delta_y) for num_peds_delta_t in 
                                          number_pass_line[line_ids.index(line_id)]])
            if not np.allclose(ref_specific_flow, out_specific_flow, atol=abs_tolerance):
                success = False
                logging.critical(f'line id <{line_id}> wrong specific flow value')
            else:
                logging.info(f'line id <{line_id}> correct specific flow value')
            
            # check density value
            if accepted_range[0] == accepted_range[1]:
                ref_density = np.array([num_peds_delta_t / (delta_t_seconds * delta_y) / accepted_range[0] for num_peds_delta_t in 
                                        number_pass_line[line_ids.index(line_id)]])
                condition = not np.allclose(ref_density, out_rho, atol=abs_tolerance)
            else:
                ref_density_0 = [num_peds_delta_t / (delta_t_seconds * delta_y) / accepted_range[1] for num_peds_delta_t in 
                                 number_pass_line[line_ids.index(line_id)]]
                ref_density_1 = [num_peds_delta_t / (delta_t_seconds * delta_y) / accepted_range[0] for num_peds_delta_t in 
                                 number_pass_line[line_ids.index(line_id)]]
                rho_valid = True
                val_idx = 0
                for rho_val in out_rho:
                    rho_valid = (ref_density_0[val_idx] <= rho_val and ref_density_1[val_idx] >= rho_val and rho_valid)
                    val_idx += 1
                condition = not rho_valid

            if condition:
                success = False
                logging.critical(f'line id <{line_id}> wrong density value')
            else:
                logging.info(f'line id <{line_id}> correct density value')

    return success

###### FUNCTIONS FOR METHOD G ########################

def get_num_peds_distance_per_dt(start_pos_x, dist, num_columns, peds_y, 
                                 dt_frames, dt_seconds, velocity, num_frames, x0, x1):
    number_time_intervals = get_num_time_intervals(num_frames, dt_frames)
    dist_per_dt = velocity * dt_seconds
    dist_per_frame = velocity * dt_seconds / dt_frames
    column_position_t0 = [start_pos_x - num_column * dist for num_column in range(num_columns)]
    column_position_t1 = [start_pos_x - num_column * dist + dist_per_dt for num_column in range(num_columns)]

    num_peds_per_dt = []
    distance_per_dt = []
    for i in range(number_time_intervals):
        num_columns_dt = 0
        sum_distance = 0
        for column_idx in range(num_columns):
            x_pos_t0 = column_position_t0[column_idx]
            x_pos_t1 = column_position_t1[column_idx]
            distance_1 = [
                math.isclose(x_pos_t0, x0, abs_tol=0.00001) and math.isclose(x_pos_t1, x1, abs_tol=0.00001),
                math.isclose(x_pos_t0, x0, abs_tol=0.00001) and x_pos_t1 > x0,
                x_pos_t0 < x1 and math.isclose(x_pos_t1, x1, abs_tol=0.00001),
                x_pos_t0 > x0 and x_pos_t1 < x1,
                ]
            # have to use isclose because of rounding deviations
            distance_2 = [
                x_pos_t0 < x0 and x_pos_t1 > x1 
                and not math.isclose(x_pos_t0, x0, abs_tol=0.00001) 
                and not math.isclose(x_pos_t1, x1, abs_tol=0.00001),
                x_pos_t0 < x0 and x_pos_t1 > x0
                and not math.isclose(x_pos_t1, x0, abs_tol=0.00001),
                x_pos_t0 < x1 and x_pos_t1 > x1
                and not math.isclose(x_pos_t0, x1, abs_tol=0.00001)
                ]
            if True in distance_1:
                num_columns_dt += 1
                sum_distance += (x_pos_t1 - x_pos_t0) * peds_y
            elif True in distance_2:
                # have to use elif because of rounding deviations
                if distance_2[0]:
                    num_columns_dt += 1
                    first_x = x_pos_t0
                    while first_x < x0:
                        first_x += dist_per_frame
                    last_x = x_pos_t1
                    while last_x > x1:
                        last_x -= dist_per_frame
                    sum_distance += (last_x - first_x) * peds_y
                elif distance_2[1]:
                    first_x = x_pos_t0
                    while first_x < x0:
                        first_x += dist_per_frame
                    if not math.isclose(first_x, x_pos_t1, abs_tol=0.00001):
                        num_columns_dt += 1
                        sum_distance += (x_pos_t1 - first_x) * peds_y
                else:
                    num_columns_dt += 1
                    last_x = x_pos_t0
                    while last_x < x1:
                        last_x += dist_per_frame
                    sum_distance += (last_x - x_pos_t0) * peds_y
        column_position_t0 = column_position_t1
        column_position_t1 = [pos + dist_per_dt for pos in column_position_t0]
        num_peds_per_dt.append(num_columns_dt * peds_y)
        distance_per_dt.append(sum_distance)

    return [num_peds_per_dt, distance_per_dt]

def get_num_peds_per_dx(start_pos_x, dist, num_columns, peds_y, 
                        delta_t_frames, delta_t_seconds, velocity, num_frames, x0, dx, num_poly):
    # x0 has to be lower boundary of x-range of MA and MA has to be not rotated

    number_time_intervals = get_num_time_intervals(num_frames, delta_t_frames)
    fps = delta_t_frames/delta_t_seconds
    dist_per_delta_t = velocity * delta_t_frames / fps

    num_peds_per_dx = []
    distance_per_dx = []
    for polygon_idx in range(num_poly):
        column_position_t0 = [start_pos_x - num_column * dist for num_column in range(num_columns)]
        column_position_t1 = [start_pos_x - num_column * dist + dist_per_delta_t for num_column in range(num_columns)]
        x0_poly = x0 + polygon_idx * dx
        x1_poly = x0 + (polygon_idx + 1) * dx

        num_columns_delta_t = []
        for i in range(number_time_intervals):
            current_num_columns = 0
            for column_idx in range(num_columns):
                x_pos_t0 = column_position_t0[column_idx]
                x_pos_t1 = column_position_t1[column_idx]

                if x_pos_t0 <= x0_poly and x_pos_t1 >= x1_poly:
                    current_num_columns += 1

            column_position_t0 = column_position_t1
            column_position_t1 = [pos + dist_per_delta_t for pos in column_position_t0]
            num_columns_delta_t.append(current_num_columns * peds_y)
        num_peds_per_dx.append(num_columns_delta_t)

    return num_peds_per_dx

def runtest_method_G(trajfile,
                     dt_frames, dt_seconds, delta_t_seconds, num_frames,
                     delta_x, n_polygon, dx, 
                     real_velocity, 
                     fps,
                     number_pass_cut_area, number_pass_area, distances_per_dt):
    success = True
    abs_tolerance = 0.0001
    delta_t_frames = delta_t_seconds * fps
    # values are rounded differently in output files -> +- 0.0001 as tolerance
    general_output_path = os.path.join('./Output', 'Fundamental_Diagram', 'Method_G')

    logging.info("===== Method G =========================")

    out_v_fname = os.path.join(general_output_path, f'v_{trajfile}_id_1.dat')
    out_rho_fname = os.path.join(general_output_path, f'rho_{trajfile}_id_1.dat')
    out_rho_flow_v_fname = os.path.join(general_output_path, f'rho_flow_v_{trajfile}_id_1.dat')

    if not os.path.exists(out_v_fname) or not os.path.exists(out_rho_fname) or not os.path.exists(out_rho_flow_v_fname):
        logging.critical("jpsreport did not output results correctly.")
        exit(FAILURE)

    # check if files were output correctly

    out_v = np.genfromtxt(out_v_fname)
    out_rho = np.genfromtxt(out_rho_fname)
    out_rho_flow_v = np.genfromtxt(out_rho_flow_v_fname)

    #### CHECK DX (FIXED PLACE) ###############

    logging.info('---- check dx values ----')

    number_time_intervals = get_num_time_intervals(num_frames, delta_t_frames)

    if not (out_v.shape == (n_polygon, number_time_intervals) and number_time_intervals > 1) and \
        not (out_v.shape == (n_polygon,) and number_time_intervals == 1):
        # should be one time interval in these scenarios
        # else the array has to be handled differently below
        success = False
        logging.critical('wrong number of velocity values for dx (fixed place)')
    elif not (out_rho.shape == (n_polygon, number_time_intervals) and number_time_intervals > 1) and \
        not (out_rho.shape == (n_polygon,) and number_time_intervals == 1):
        # should be one time interval in these scenarios
        # else the array has to be handled differently below
        success = False
        logging.critical('wrong number of density values for dx (fixed place)')
    else:
        logging.info('correct number of values for dx (fixed place)')
        
        v_range_dx = get_velocity_range(dx, real_velocity, fps, abs_tolerance)
        # the same reason for v_range as in method_F
        # however: v_range is not needed for velocity calculation over dt
        # as pedestrians do not have to pass the whole area to be included in
        # velocity calculation -> the velocity is calculated for each dt
        # (also distances are not discrete, but frames are)

        # iterate though all polygons
        idx_polygon = 0
        for v_values in out_v:
            logging.info(f'---- cut polygon {idx_polygon} ----')

            if number_time_intervals == 1:
                v_values = [v_values]

            # iterate though all time intervals
            for idx_time_interval in range(number_time_intervals):
                # check velocity value
                if number_pass_cut_area[idx_polygon][idx_time_interval] == 0:
                    condition = not np.isnan(v_values[idx_time_interval])
                elif v_range_dx[0] == v_range_dx[1]:
                    condition = not math.isclose(v_values[idx_time_interval], v_range_dx[0], abs_tol=abs_tolerance)
                else:
                    condition = v_range_dx[0] > v_values[idx_time_interval] or v_values[idx_time_interval] > v_range_dx[1]
                if condition:
                    success = False
                    logging.critical(f'velocity value {v_values[idx_time_interval]} for dx (fixed place) is not in accepted '
                                     f'range around {real_velocity} m/s')
                else:
                    logging.info('correct velocity value for dx (fixed place)')

                # check density value
                if number_time_intervals == 1:
                    rho_value = out_rho[idx_polygon]
                else:
                    rho_value = out_rho[idx_polygon][idx_time_interval]
                density_formula = (number_pass_cut_area[idx_polygon][idx_time_interval] * dx / real_velocity) / (delta_t_seconds * dx)
                # velocity formula -> rewriting the formula and inserting the real velocity to get the sum of the time
                # -> might not represent the actual density, but the density according to the formula used in method G 
                # -> e.g. real density might be 10, but 9.92555 here (which is also what jpsreport should output)
                # reason: last frame(s) in which no pedestrian passes the whole area (or a whole polygon)
                # to get the actual density one would have to subtract these frame(s) or this part of a frame
                # from delta_t_seconds (after converting it to seconds)
                # this is the reason why the formula is used to check results

                rho_range_dx = get_density_range(dx, real_velocity, fps, abs_tolerance, number_pass_cut_area[idx_polygon][idx_time_interval], delta_t_seconds)
                if rho_range_dx[0] == rho_range_dx[1]:
                    condition = not math.isclose(rho_value, rho_range_dx[0], abs_tol=abs_tolerance)
                else:
                    condition = rho_range_dx[0] > rho_value or rho_value > rho_range_dx[1]
                if condition:
                    success = False
                    logging.critical(f'density value {rho_value} for dx (fixed place) is not in accepted '
                                        f'range around {density_formula} 1/m')
                else:
                    logging.info('correct density value for dx (fixed place)')

            idx_polygon += 1

    #### CHECK DT (FIXED TIME) ###############
    
    logging.info('---- check dt values ----')
    number_time_intervals = get_num_time_intervals(num_frames - (num_frames % delta_t_frames) + 1, dt_frames)
    if out_rho_flow_v.shape[0] != number_time_intervals:
        # check number of time intervals
        success = False
        logging.critical('wrong number of time intervals for dt (fixed time)')
    else:
        logging.info('correct number of time intervals for dt (fixed time)')

        # check all velocity, density and flow values
        dt_idx = 0
        reference_v = []
        reference_flow = []
        for num_peds_dt in number_pass_area:
            if num_peds_dt != 0:
                reference_v.append(distances_per_dt[dt_idx] / (dt_seconds * num_peds_dt))
            else:
                reference_v.append(np.nan)
            reference_flow.append(distances_per_dt[dt_idx] / (dt_seconds * delta_x))
            dt_idx += 1
        reference_v = np.array(reference_v[:number_time_intervals])
        reference_flow = np.array(reference_flow[:number_time_intervals])
        reference_rho = np.array([num_peds_dt / delta_x for num_peds_dt in number_pass_area][:number_time_intervals])

        if not np.allclose(reference_v, out_rho_flow_v[:, 0], atol=abs_tolerance, equal_nan=True):
            success = False
            logging.critical(f'wrong velocity values for dt (fixed time); real velocity {real_velocity} m/s')
        else:
            logging.info('correct velocity values for dt (fixed time)')

        if not np.allclose(reference_rho, out_rho_flow_v[:, 1], atol=abs_tolerance):
            success = False
            logging.critical('wrong density values for dt (fixed time)')
        else:
            logging.info('correct density values for dt (fixed time)')

        if not np.allclose(reference_flow, out_rho_flow_v[:, 2], atol=abs_tolerance):
            success = False
            logging.critical('wrong flow values for dt (fixed time)')
        else:
            logging.info('correct flow values for dt (fixed time)')

    return success

###### FUNCTIONS FOR METHOD H ########################

def runtest_method_H(trajfile,
                     delta_t_frames, delta_t_seconds, num_frames,
                     delta_x, 
                     real_velocity, 
                     fps,
                     accum_peds_delta_t):
    success = True
    abs_tolerance = 0.0001
    # values are rounded differently in output files -> +- 0.0001 as tolerance
    general_output_path = os.path.join('./Output', 'Fundamental_Diagram', 'Method_H')

    logging.info("===== Method H =========================")

    out_fname = os.path.join(general_output_path, f'flow_rho_v_{trajfile}_id_1.dat')
    if not os.path.exists(out_fname):
        # check if file was output correctly
        logging.critical("jpsreport did not output results correctly.")
        exit(FAILURE)

    out_data = np.genfromtxt(out_fname)

    number_time_intervals = get_num_time_intervals(num_frames, delta_t_frames)

    if not (out_data.shape == (number_time_intervals, 3) and number_time_intervals > 1) and \
            not (out_data.shape == (3,) and number_time_intervals == 1):
        success = False
        logging.critical('wrong number of time intervals')
    else:
        logging.info('correct number of time intervals')

        if number_time_intervals == 1:
            ref_rho = accum_peds_delta_t[0] / delta_t_frames / delta_x
            # accum_peds_delta_t is the sum of the pedestrians that are on the MA at a frame
            # for each frame -> is "sum of time"
            # delta_t_frames does not have to be changed to seconds, as the ratio
            # between accum_peds_delta_t and delta_t_frames matters
            # -> accum_peds_delta_t / delta_t_frames is the average number of peds
            # in MA at one frame

            ref_flow = ref_rho * real_velocity

            condition_flow = not math.isclose(ref_flow, out_data[0], abs_tol=abs_tolerance)
            condition_rho = not math.isclose(ref_rho, out_data[1], abs_tol=abs_tolerance)
            condition_v = not math.isclose(real_velocity, out_data[2], abs_tol=abs_tolerance)
        else:
            ref_rho = np.array([num_peds / delta_t_frames / delta_x for num_peds in accum_peds_delta_t])
            ref_flow = ref_rho * real_velocity

            condition_flow = not np.allclose(ref_flow, out_data[:, 0], atol=abs_tolerance)
            condition_rho = not np.allclose(ref_rho, out_data[:, 1], atol=abs_tolerance)
            condition_v = not np.allclose(real_velocity, out_data[:, 2], atol=abs_tolerance)

        if condition_flow:
            success = False
            logging.critical("wrong flow values")
        else:
            logging.info("correct flow values")

        if condition_rho:
            success = False
            logging.critical("wrong density values")
        else:
            logging.info("correct density values")

        if condition_v:
            success = False
            logging.critical("wrong velocity values")
        else:
            logging.info("correct velocity values")

    return success