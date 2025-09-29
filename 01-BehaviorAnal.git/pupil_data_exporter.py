# -*- coding: utf-8 -*-
# Copyright (c) 2022. Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS
# @Time     : 2022/11/1 19:40 
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: https://github.com/pupil-labs/pupil/tree/master/pupil_src/shared_modules
# https://gist.github.com/papr/f7d5fbf7809578e26e1cbb77fdb9ad06
# https://github.com/pupil-labs/pupil-community
# https://github.com/PierreYvesH/Pupil_labs_batch_exporter/blob/e2dfefb19752a30904f6b66d67b4c916321b33e0/from_pupil/extract_diameter.py
# @FileName : pupil_data_exporter.py
# @Software : Python3.6; PyCharm; Windows10 / Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-79-generic x86_64)
# @Hardware : Intel Core i7-4712MQ; NVIDIA GeForce 840M / 2*X640-G30(XEON 6258R 2.7G); 3*NVIDIA GeForce RTX3090
# @Version  : V1.1 - ZL.Z：2023/5/3
# 		      处理某些文件不存在报错问题
#             V1.0 - ZL.Z：2022/11/1 - 2022/11/3, 2022/11/6
# 		      First version.
# @License  : None
# @Brief    : 提取并保存Pupil数据，和Pupil Player导出的数据几乎一致

import logging
import os
import numpy as np
import msgpack
import pandas as pd
import json
from scipy.signal import fftconvolve
from scipy.spatial.distance import pdist
import typing
from typing import NamedTuple
from collections import deque
from functools import partial
from pathos.pools import ProcessPool as Pool
import argparse

logger = logging.getLogger(__name__)


class _BaseDetectionExporter:
    """各类Pupil的检测数据提取基类"""
    def __init__(self, topic: str = "pupil", min_conf_thr: float = 0.0):
        """
        :param topic: 待提取的主题名
        :param min_conf_thr: 最小置信阈值
                            Data with a lower confidence than this threshold is not considered during the detection.
        """
        self.topic = topic
        self.csv_export_filename = self.topic + '.csv'
        if min_conf_thr is None:
            min_conf_thr = 0.0
        self.min_conf_thr = min_conf_thr

    @classmethod
    def dict_export(cls, raw_value: dict, world_index: int) -> dict:
        """
        获取字典值
        :param raw_value: 解包后的pupil datum数据
        :param world_index: 世界视频帧索引
        :return: 待导出的字典值
        """
        pass

    def process_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                           overwrite=False, n_jobs=None) -> typing.Union[pd.DataFrame, list]:
        """
        处理单个或多个瞳孔测量记录（多个时，可选并行或非并行）
        :param recordings: 处理单/多个记录，可以是包含瞳孔记录文件的单个文件夹，
                           也可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 处理后的所有数据 pd.DataFrame类型（单个瞳孔测量记录时）或list类型（多个瞳孔测量记录时）
        """
        if type(recordings) == list:
            res = self.process_multi_recordings(recordings, csv_out_dir, overwrite, n_jobs)
        else:
            rec_l = []
            for recording in os.listdir(recordings):
                if os.path.isfile(os.path.join(recordings, recording)):
                    rec_l = []
                    break
                rec_l.append(os.path.join(recordings, recording))
            if rec_l:
                res = self.process_multi_recordings(recordings, csv_out_dir, overwrite, n_jobs)
            else:
                res = self.process_single_recording(recordings, csv_out_dir, overwrite)
        return res

    def process_single_recording(self, recording: str, csv_out_dir: str = '',
                                 overwrite=False) -> typing.Optional[pd.DataFrame]:
        """
        处理单个瞳孔测量记录
        :param recording: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径或为带文件名的路径）
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :return: 所有数据 pd.DataFrame类型
        """
        print(f'---------- Processing {recording} ... ----------')
        if csv_out_dir == '':
            export_f = os.path.join(recording, 'exports', os.path.basename(recording), self.csv_export_filename)
        elif csv_out_dir.endswith('.csv'):
            export_f = csv_out_dir
        else:
            export_f = os.path.join(csv_out_dir, os.path.basename(recording.rstrip('/')), self.csv_export_filename)
        if os.path.exists(export_f):
            if not overwrite:
                print("{} exists already! Not overwriting, skipping.".format(export_f))
                return None
            else:
                logger.warning("{} exists already! Overwriting.".format(export_f))
        sort_by = self.topic + "_timestamp"
        try:
            exporter_data = pd.DataFrame(self.load_and_yield_data(recording))
        except FileNotFoundError as e:
            print(f"文件缺失/格式有误，略过：{e}")
            return None
        if exporter_data.empty:
            raise ValueError(f"No {self.topic} data found, {self.topic}.pldata file may be empty, please check!")
        exporter_data.sort_values(by=sort_by, ignore_index=True, inplace=True)
        if not os.path.exists(os.path.dirname(export_f)):
            os.makedirs(os.path.dirname(export_f))
        exporter_data.to_csv(export_f, index=False, encoding='utf-8-sig')
        return exporter_data

    def process_multi_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                                 overwrite=False, n_jobs=None) -> list:
        """
        处理多个瞳孔测量记录（并行或非并行）
        :param recordings: 处理多个记录，可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 并行处理后的所有数据 list类型
        """
        assert os.path.isdir(csv_out_dir) or csv_out_dir == '', "csv_out_dir仅接受已存在路径或空值''输入"
        assert (n_jobs is None) or (type(n_jobs) is int and n_jobs > 0) or (
                n_jobs == -1), 'n_jobs仅接受-1/正整数/None类型输入'
        if n_jobs == -1:
            n_jobs = None
        rec_l = []
        if type(recordings) == list:
            rec_l = recordings
        else:
            for recording in os.listdir(recordings):
                if os.path.isfile(os.path.join(recordings, recording)):
                    raise ValueError(u'输入文件夹内存在其他数据文件，请重新选择！')
                rec_l.append(os.path.join(recordings, recording))
        if n_jobs == 1:
            res = []
            for recording in rec_l:
                res.append(self.process_single_recording(recording, csv_out_dir, overwrite))
        else:
            parallel_process = partial(self.process_single_recording, csv_out_dir=csv_out_dir, overwrite=overwrite)
            with Pool(n_jobs) as pool:
                res = pool.map(parallel_process, rec_l)
        return res

    def load_and_yield_data(self, directory: str) -> typing.Generator:
        """
        Load and extract pupil diameter data
        :param directory: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :return: <class 'generator'> 所提取到的所有数据的生成器
        """
        data_ts = np.load(os.path.join(directory, self.topic + "_timestamps.npy"))
        with open(os.path.join(directory, "info.player.json"), 'r') as f:
            info_data = json.load(f)
        try:
            world_ts = np.load(os.path.join(directory, "world_timestamps.npy"))
        except FileNotFoundError:
            logger.warning("world_timestamps.npy not exist already, "
                           "using the fake timestamps with 30Hz sampling rate of the world camera!")
            world_ts = np.arange(info_data['start_time_synced_s'],
                                 info_data['start_time_synced_s'] + info_data['duration_s'], 1/30.0)
        export_world_idc = find_closest(world_ts, data_ts)  # 获取最近的世界视频帧
        msgpack_file = os.path.join(directory, self.topic + ".pldata")
        with open(msgpack_file, "rb") as fh:
            unpacker = msgpack.Unpacker(fh, raw=False, use_list=False)
            for ts, (topic, payload), idx in zip(data_ts, unpacker, export_world_idc):
                datum = msgpack.unpackb(payload, raw=False, use_list=False)
                if datum["confidence"] < self.min_conf_thr or (ts < world_ts[0]) or (ts > world_ts[-1]):
                    continue
                dict_row = self.dict_export(datum, idx)
                yield dict_row


class PupilDetectionExporter(_BaseDetectionExporter):
    """瞳孔检测数据提取"""
    def __init__(self, topic: str = "pupil", **kwargs):
        """
        :param topic: 待提取的主题名
        :param **kwargs: 父类__init__参数
        """
        super().__init__(topic=topic, **kwargs)

    @classmethod
    def dict_export(cls, raw_value: dict, world_index: int) -> dict:
        """
        获取字典值
        :param raw_value: 解包后的pupil datum数据
        :param world_index: 世界视频帧索引
        :return: 待导出的字典值
        """
        # 2d data
        pupil_timestamp = raw_value["timestamp"]
        eye_id = raw_value["id"]
        confidence = raw_value["confidence"]
        norm_pos_x = raw_value["norm_pos"][0]
        norm_pos_y = raw_value["norm_pos"][1]
        diameter = raw_value["diameter"]
        method = raw_value["method"]
        # ellipse data
        try:
            ellipse_center = raw_value["ellipse"]["center"]
            ellipse_axis = raw_value["ellipse"]["axes"]
            ellipse_angle = raw_value["ellipse"]["angle"]
        except KeyError:
            ellipse_center = [None, None]
            ellipse_axis = [None, None]
            ellipse_angle = None
        # 3d data
        try:
            diameter_3d = raw_value["diameter_3d"]
            model_confidence = raw_value["model_confidence"]
            sphere_center = raw_value["sphere"]["center"]
            sphere_radius = raw_value["sphere"]["radius"]
            circle_3d_center = raw_value["circle_3d"]["center"]
            circle_3d_normal = raw_value["circle_3d"]["normal"]
            circle_3d_radius = raw_value["circle_3d"]["radius"]
            theta = raw_value["theta"]
            phi = raw_value["phi"]
            projected_sphere_center = raw_value["projected_sphere"]["center"]
            projected_sphere_axis = raw_value["projected_sphere"]["axes"]
            projected_sphere_angle = raw_value["projected_sphere"]["angle"]
        except KeyError:
            diameter_3d = None
            model_confidence = None
            sphere_center = [None, None, None]
            sphere_radius = None
            circle_3d_center = [None, None, None]
            circle_3d_normal = [None, None, None]
            circle_3d_radius = None
            theta = None
            phi = None
            projected_sphere_center = [None, None]
            projected_sphere_axis = [None, None]
            projected_sphere_angle = None
        return {
            # 2d data
            "pupil_timestamp": pupil_timestamp,
            "world_index": world_index,
            "eye_id": eye_id,
            "confidence": confidence,
            "norm_pos_x": norm_pos_x,
            "norm_pos_y": norm_pos_y,
            "diameter": diameter,
            "method": method,
            # ellipse data
            "ellipse_center_x": ellipse_center[0],
            "ellipse_center_y": ellipse_center[1],
            "ellipse_axis_a": ellipse_axis[0],
            "ellipse_axis_b": ellipse_axis[1],
            "ellipse_angle": ellipse_angle,
            # 3d data
            "diameter_3d": diameter_3d,
            "model_confidence": model_confidence,
            "sphere_center_x": sphere_center[0],
            "sphere_center_y": sphere_center[1],
            "sphere_center_z": sphere_center[2],
            "sphere_radius": sphere_radius,
            "circle_3d_center_x": circle_3d_center[0],
            "circle_3d_center_y": circle_3d_center[1],
            "circle_3d_center_z": circle_3d_center[2],
            "circle_3d_normal_x": circle_3d_normal[0],
            "circle_3d_normal_y": circle_3d_normal[1],
            "circle_3d_normal_z": circle_3d_normal[2],
            "circle_3d_radius": circle_3d_radius,
            "theta": theta,
            "phi": phi,
            "projected_sphere_center_x": projected_sphere_center[0],
            "projected_sphere_center_y": projected_sphere_center[1],
            "projected_sphere_axis_a": projected_sphere_axis[0],
            "projected_sphere_axis_b": projected_sphere_axis[1],
            "projected_sphere_angle": projected_sphere_angle,
        }


class GazeDetectionExporter(_BaseDetectionExporter):
    """凝视检测数据提取"""
    def __init__(self, topic: str = "gaze", **kwargs):
        """
        :param topic: 待提取的主题名
        :param **kwargs: 父类__init__参数
        """
        super().__init__(topic=topic, **kwargs)

    @classmethod
    def dict_export(cls, raw_value, world_index: int) -> dict:
        """
        获取字典值
        :param raw_value: 解包后的pupil datum数据
        :param world_index: 世界视频帧索引
        :return: 待导出的字典值
        """
        gaze_timestamp = raw_value["timestamp"]
        confidence = raw_value["confidence"]
        norm_pos = raw_value["norm_pos"]
        base_data = None
        gaze_points_3d = [None, None, None]
        eye_centers0_3d = [None, None, None]
        eye_centers1_3d = [None, None, None]
        gaze_normals0_3d = [None, None, None]
        gaze_normals1_3d = [None, None, None]
        if raw_value.get("base_data", None) is not None:
            base_data = " ".join("{}-{}".format(b["timestamp"], b["id"]) for b in raw_value["base_data"])
        # add 3d data if avaiblable
        if raw_value.get("gaze_point_3d", None) is not None:
            gaze_points_3d = raw_value["gaze_point_3d"]
            # binocular
            if raw_value.get("eye_centers_3d", None) is not None:
                eye_centers_3d = raw_value["eye_centers_3d"]
                gaze_normals_3d = raw_value["gaze_normals_3d"]
                eye_centers0_3d = (
                    eye_centers_3d.get("0", None)
                    or eye_centers_3d.get(0, None)  # backwards compatibility
                    or [None, None, None]
                )
                eye_centers1_3d = (
                    eye_centers_3d.get("1", None)
                    or eye_centers_3d.get(1, None)  # backwards compatibility
                    or [None, None, None]
                )
                gaze_normals0_3d = (
                    gaze_normals_3d.get("0", None)
                    or gaze_normals_3d.get(0, None)  # backwards compatibility
                    or [None, None, None]
                )
                gaze_normals1_3d = (
                    gaze_normals_3d.get("1", None)
                    or gaze_normals_3d.get(1, None)  # backwards compatibility
                    or [None, None, None]
                )
            # monocular
            elif raw_value.get("eye_center_3d", None) is not None:
                try:
                    eye_id = raw_value["base_data"][0]["id"]
                except (KeyError, IndexError):
                    logger.warning(
                        f"Unexpected raw base_data for monocular gaze!"
                        f" Data: {raw_value.get('base_data', None)}"
                    )
                else:
                    if str(eye_id) == "0":
                        eye_centers0_3d = raw_value["eye_center_3d"]
                        gaze_normals0_3d = raw_value["gaze_normal_3d"]
                    elif str(eye_id) == "1":
                        eye_centers1_3d = raw_value["eye_center_3d"]
                        gaze_normals1_3d = raw_value["gaze_normal_3d"]
        return {
            "gaze_timestamp": gaze_timestamp,
            "world_index": world_index,
            "confidence": confidence,
            "norm_pos_x": norm_pos[0],
            "norm_pos_y": norm_pos[1],
            "base_data": base_data,
            "gaze_point_3d_x": gaze_points_3d[0],
            "gaze_point_3d_y": gaze_points_3d[1],
            "gaze_point_3d_z": gaze_points_3d[2],
            "eye_center0_3d_x": eye_centers0_3d[0],
            "eye_center0_3d_y": eye_centers0_3d[1],
            "eye_center0_3d_z": eye_centers0_3d[2],
            "gaze_normal0_x": gaze_normals0_3d[0],
            "gaze_normal0_y": gaze_normals0_3d[1],
            "gaze_normal0_z": gaze_normals0_3d[2],
            "eye_center1_3d_x": eye_centers1_3d[0],
            "eye_center1_3d_y": eye_centers1_3d[1],
            "eye_center1_3d_z": eye_centers1_3d[2],
            "gaze_normal1_x": gaze_normals1_3d[0],
            "gaze_normal1_y": gaze_normals1_3d[1],
            "gaze_normal1_z": gaze_normals1_3d[2],
        }


class _FixationDatum(NamedTuple):
    """Algorithm output"""
    id: int
    start_timestamp: float
    duration_ms: float
    start_frame_index: int
    end_frame_index: int
    norm_pos_x_avg: float
    norm_pos_y_avg: float
    dispersion_deg: float
    confidence_avg: float
    method: str
    gaze_point_3d_x_avg: float
    gaze_point_3d_y_avg: float
    gaze_point_3d_z_avg: float
    base_data: str


class _FixationResultFactory(object):

    def __init__(self, directory: str):
        """
        工厂模式初始化
        :param directory: 处理单个记录，包含瞳孔记录文件的单个文件夹
        """
        self._id_counter = 0
        with open(os.path.join(directory, "info.player.json"), 'r') as f:
            info_data = json.load(f)
        try:
            self.world_ts = np.load(os.path.join(directory, "world_timestamps.npy"))
        except FileNotFoundError:
            logger.warning("world_timestamps.npy not exist already, "
                           "using the fake timestamps with 30Hz sampling rate of the world camera!")
            self.world_ts = np.arange(info_data['start_time_synced_s'],
                                      info_data['start_time_synced_s'] + info_data['duration_s'], 1/30.0)

    def from_data(self, dispersion_result: float,
                  data: typing.Union[typing.Deque, list]) -> _FixationDatum:
        """
        计算获取注视检测数据
        :param dispersion_result: Dispersion, in degrees
        :param data: 每个作为注视时间片段的凝视检测数据
        :return: _FixationDatum类型的注视检测数据
        """
        df = pd.DataFrame(data)
        start_ts, stop_ts = df.timestamp.iloc[[0, -1]]
        duration_ms = (stop_ts - start_ts) * 1000
        direction_3d = df[list("xyz")].mean(axis=0)
        norm_2d = df.norm_pos.mean(axis=0)
        base_data = []
        for i_base in df.base_data:
            for j_base in i_base:
                base_data.append(j_base["timestamp"])
        # correlate world indices
        start, end = np.searchsorted(self.world_ts, [base_data[0], base_data[-1]])
        end = min(end, len(self.world_ts) - 1)  # fix `list index out of range` error
        fixation = _FixationDatum(
            id=self._id_counter,
            start_timestamp=start_ts,
            duration_ms=duration_ms,
            start_frame_index=int(start),
            end_frame_index=int(end),
            dispersion_deg=dispersion_result,
            norm_pos_x_avg=norm_2d[0],
            norm_pos_y_avg=norm_2d[1],
            confidence_avg=df.confidence.mean(axis=0),
            method='3d gaze' if df.gp_3d_flag.all() else '2d gaze',
            gaze_point_3d_x_avg=direction_3d[0],
            gaze_point_3d_y_avg=direction_3d[1],
            gaze_point_3d_z_avg=direction_3d[2],
            base_data=" ".join(str(b) for b in base_data),
        )
        self._id_counter += 1
        return fixation


class FixationDetectionExporter:
    """注视检测数据提取"""
    def __init__(self, max_dispersion_deg: float = 1.50, min_duration_ms: int = 100,
                 max_duration_ms: int = 400, min_conf_thr: float = 0.6):
        """
        阈值参考：max_dispersion_deg=1.5(Pupil Player 默认值)；min_duration_ms=100；max_duration_ms=400；
        min_conf_thr=0.6(Pupil Player 默认值)
        1. Salvucci, D. D., & Goldberg, J. H. (2000). Identifying fixations and saccades in eye-tracking protocols.
        Paper presented at the Proceedings of the symposium on Eye tracking research & applications - ETRA '00.
        2. Manor, B. R., & Gordon, E. (2003). Defining the temporal threshold for ocular fixation in free-viewing
        visuocognitive tasks. Journal of Neuroscience Methods, 128(1-2), 85-93. doi:10.1016/s0165-0270(03)00151-1
        3. Blignaut, P. (2009). Fixation identification: the optimum threshold for a dispersion algorithm. Atten
        Percept Psychophys, 71(4), 881-895. doi:10.3758/APP.71.4.881
        4. Orquin, J. L., & Holmqvist, K. (2018). Threats to the validity of eye-movement research in psychology.
        Behav. Res. Methods, 50(4), 1645-1656. doi:10.3758/s13428-017-0998-z
        5. https://www.researchgate.net/post/Fixation_duration_threshold_of_cognitive_processing
        :param max_dispersion_deg: Maximum distance between all gaze locations during a fixation.
        :param min_duration_ms: The minimum duration in which the dispersion threshold must not be exceeded.
        :param max_duration_ms: The maximum duration in which the dispersion threshold must not be exceeded.
        :param min_conf_thr: Data with a lower confidence than this threshold is not
                             considered during fixation detection.
        """
        self.csv_export_filename = "fixations.csv"
        self.max_dispersion_deg = max_dispersion_deg
        self.min_duration_ms = min_duration_ms
        self.max_duration_ms = max_duration_ms
        if min_conf_thr is None:
            min_conf_thr = 0.6
        self.min_conf_thr = min_conf_thr

    @classmethod
    def dict_export(cls, raw_value: dict, intrinsics: dict) -> dict:
        """
        获取字典值
        :param raw_value: 解包后的pupil datum数据
        :param intrinsics: intrinsics数据
        :return: 待导出的字典值
        """
        timestamp = raw_value["timestamp"]
        confidence = raw_value["confidence"]
        base_data = None
        # add 3d data if avaiblable
        if raw_value.get("gaze_point_3d", None) is not None:
            gaze_points_3d = raw_value["gaze_point_3d"]
            gp_3d_flag = True
        else:
            gaze_points_3d = normalized2d_to_direction3d(raw_value["norm_pos"], intrinsics)
            gp_3d_flag = False
        if raw_value.get("base_data", None) is not None:
            base_data = raw_value["base_data"]
        return {
            "timestamp": timestamp,
            "confidence": confidence,
            "gp_3d_flag": gp_3d_flag,
            'norm_pos': np.array(raw_value["norm_pos"]),
            "x": gaze_points_3d[0],
            "y": gaze_points_3d[1],
            "z": gaze_points_3d[2],
            "base_data": base_data,
        }

    def process_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                           overwrite=False, n_jobs=None) -> typing.Union[pd.DataFrame, list]:
        """
        处理单个或多个瞳孔测量记录（多个时，可选并行或非并行）
        :param recordings: 处理单/多个记录，可以是包含瞳孔记录文件的单个文件夹，
                           也可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 处理后的所有数据 pd.DataFrame类型（单个瞳孔测量记录时）或list类型（多个瞳孔测量记录时）
        """
        if type(recordings) == list:
            res = self.process_multi_recordings(recordings, csv_out_dir, overwrite, n_jobs)
        else:
            rec_l = []
            for recording in os.listdir(recordings):
                if os.path.isfile(os.path.join(recordings, recording)):
                    rec_l = []
                    break
                rec_l.append(os.path.join(recordings, recording))
            if rec_l:
                res = self.process_multi_recordings(recordings, csv_out_dir, overwrite, n_jobs)
            else:
                res = self.process_single_recording(recordings, csv_out_dir, overwrite)
        return res

    def process_single_recording(self, recording: str, csv_out_dir: str = '', overwrite=False) -> typing.Optional[pd.DataFrame]:
        """
        处理单个瞳孔测量记录
        :param recording: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径或为带文件名的路径）
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :return: 所有数据 pd.DataFrame类型
        """
        if csv_out_dir == '':
            export_f = os.path.join(recording, 'exports', os.path.basename(recording), self.csv_export_filename)
        elif csv_out_dir.endswith('.csv'):
            export_f = csv_out_dir
        else:
            export_f = os.path.join(csv_out_dir, os.path.basename(recording.rstrip('/')), self.csv_export_filename)
        if os.path.exists(export_f):
            if not overwrite:
                print("{} exists already! Not overwriting, skipping.".format(export_f))
                return None
            else:
                logger.warning("{} exists already! Overwriting.".format(export_f))
        intrinsics = self.read_intrinsics(recording)
        try:
            gaze_data = pd.DataFrame(self.load_and_yield_data(recording, intrinsics))
            if gaze_data.empty:
                raise ValueError("No gaze data found, gaze.pldata file may be empty, please check!")
            gaze_data.sort_values(by="timestamp", ignore_index=True, inplace=True)
            fixation_data = self.fixation_detection(recording, gaze_data)
            if not os.path.exists(os.path.dirname(export_f)):
                os.makedirs(os.path.dirname(export_f))
            fixation_data.to_csv(export_f, index=False, encoding='utf-8-sig')
            fixation_report = pd.DataFrame({'fixation_classifier': ['Dispersion_Duration'],
                                            'max_dispersion_deg': [self.max_dispersion_deg],
                                            'min_duration_ms': [self.min_duration_ms],
                                            'max_duration_ms': [self.max_duration_ms],
                                            'min_conf_thr': [self.min_conf_thr],
                                            'fixation_count': [len(fixation_data)]})
            fixation_report.to_csv(os.path.join(os.path.dirname(export_f), 'fixation_report.csv'),
                                   index=False, encoding='utf-8-sig')
            return fixation_data
        except FileNotFoundError as e:
            print(f"文件缺失/格式有误，略过：{e}")
            return None

    def process_multi_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                                 overwrite=False, n_jobs=None) -> list:
        """
        处理多个瞳孔测量记录（并行或非并行）
        :param recordings: 处理多个记录，可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 并行处理后的所有数据 list类型
        """
        assert os.path.isdir(csv_out_dir) or csv_out_dir == '', "csv_out_dir仅接受已存在路径或空值''输入"
        assert (n_jobs is None) or (type(n_jobs) is int and n_jobs > 0) or (
                n_jobs == -1), 'n_jobs仅接受-1/正整数/None类型输入'
        if n_jobs == -1:
            n_jobs = None
        rec_l = []
        if type(recordings) == list:
            rec_l = recordings
        else:
            for recording in os.listdir(recordings):
                if os.path.isfile(os.path.join(recordings, recording)):
                    raise ValueError(u'输入文件夹内存在其他数据文件，请重新选择！')
                rec_l.append(os.path.join(recordings, recording))
        if n_jobs == 1:
            res = []
            for recording in rec_l:
                res.append(self.process_single_recording(recording, csv_out_dir, overwrite))
        else:
            parallel_process = partial(self.process_single_recording, csv_out_dir=csv_out_dir, overwrite=overwrite)
            with Pool(n_jobs) as pool:
                res = pool.map(parallel_process, rec_l)
        return res

    def load_and_yield_data(self, directory: str, intrinsics: dict) -> typing.Generator:
        """
        Load and extract pupil diameter data
        :param directory: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :param intrinsics: intrinsics数据
        :return: <class 'generator'> 所提取到的所有数据的生成器
        """
        msgpack_file = os.path.join(directory, "gaze.pldata")
        with open(msgpack_file, "rb") as fh:
            unpacker = msgpack.Unpacker(fh, raw=False, use_list=False)
            for topic, payload in unpacker:
                datum = msgpack.unpackb(payload, raw=False, use_list=False)
                if datum["confidence"] < self.min_conf_thr:
                    continue
                dict_row = self.dict_export(datum, intrinsics)
                yield dict_row

    @classmethod
    def read_intrinsics(cls, directory: str) -> dict:
        """
        读取.intrinsics文件
        :param directory: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :return: intrinsics数据 dict类型
        """
        intrinsics_location = os.path.join(directory, "world.intrinsics")
        try:
            with open(intrinsics_location, "rb") as fh:
                intrinsics = msgpack.unpack(fh, raw=False, use_list=False)
            intrinsics = next(values for key, values in intrinsics.items() if key != "version")
        except FileNotFoundError:
            logger.warning("world.intrinsics not exist already, using the default intrinsics with "
                           "'Pupil Cam1 ID2' (1280, 720) resolution!")
            intrinsics = default_camera_intrinsics()['Pupil Cam1 ID2']['(1280, 720)']
        if intrinsics["cam_type"] not in ("radial", "fisheye"):
            raise ValueError(f"Unexpected camera model {intrinsics['cam_type']}")
        return intrinsics

    def fixation_detection(self, directory: str, gaze_data: pd.DataFrame) -> pd.DataFrame:
        """
        注视点检测
        ref: https://github.com/pupil-labs/pupil/blob/master/pupil_src/shared_modules/fixation_detector.py#L162-L253
        :param directory: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :param gaze_data: 凝视检测数据
        :return: 所有注视点检测数据 pd.DataFrame类型
        """
        fixation_result = _FixationResultFactory(directory)
        working_slice = deque()
        _, rows = zip(*gaze_data.iterrows())
        remaining_slice = deque(rows)
        fixations = []
        while remaining_slice:
            # check if working_queue contains enough data
            if len(working_slice) < 2 or \
                    (working_slice[-1]["timestamp"] - working_slice[0]["timestamp"]) < self.min_duration_ms / 1000.0:
                datum = remaining_slice.popleft()
                working_slice.append(datum)
                continue
            # min duration reached, check for fixation
            dispersion = vector_dispersion(working_slice)
            if dispersion > self.max_dispersion_deg:
                # not a fixation, move forward
                working_slice.popleft()
                continue
            left_idx = len(working_slice)
            # minimal fixation found. collect maximal data
            # to perform binary search for fixation end
            while remaining_slice:
                datum = remaining_slice[0]
                if datum["timestamp"] > working_slice[0]["timestamp"] + self.max_duration_ms / 1000.0:
                    break  # maximum data found
                working_slice.append(remaining_slice.popleft())
            # check for fixation with maximum duration
            dispersion = vector_dispersion(working_slice)
            if dispersion <= self.max_dispersion_deg:
                fixation = fixation_result.from_data(dispersion, working_slice)
                if fixation.end_frame_index:
                    fixations.append(fixation)
                working_slice.clear()  # discard old Q
                continue
            slicable = list(working_slice)  # deque does not support slicing
            right_idx = len(working_slice)
            # binary search
            while left_idx < right_idx - 1:
                middle_idx = (left_idx + right_idx) // 2
                dispersion = vector_dispersion(slicable[: middle_idx + 1])
                if dispersion <= self.max_dispersion_deg:
                    left_idx = middle_idx
                else:
                    right_idx = middle_idx
            # left_idx-1 is last valid base datum
            final_base_data = slicable[:left_idx]
            to_be_placed_back = slicable[left_idx:]
            dispersion_result = vector_dispersion(final_base_data)
            fixation = fixation_result.from_data(dispersion_result, final_base_data)
            if fixation.end_frame_index:
                fixations.append(fixation)
            working_slice.clear()  # clear queue
            remaining_slice.extendleft(reversed(to_be_placed_back))
        fixations_data = pd.DataFrame(fixations)
        return fixations_data


class _BlinkDatum(NamedTuple):
    id: int
    start_timestamp: float
    duration_s: float
    end_timestamp: float
    start_frame_index: int
    index: int
    end_frame_index: int
    confidence: float
    filter_response: str
    base_data: str


class BlinkDetectionExporter:
    """眨眼检测数据提取"""
    def __init__(self, filter_len_s: float = 0.20, onset_conf_thr: float = 0.50, offset_conf_thr: float = 0.50):
        """
        阈值参考(Pupil Player 默认值)：filter_len_s=0.20, onset_conf_thr=0.50, offset_conf_thr=0.50
        :param filter_len_s: The time window's length in which the detector tries to find confidence drops and gains.
        :param onset_conf_thr: The threshold that the filter response ('Activity') must rise above to classify
                               the onset of a blink, corresponding to a sudden drop in 2D pupil detection confidence.
        :param offset_conf_thr: The threshold that the filter response ('Activity') must fall below to classify
                                the end of a blink, corresponding to a rise in 2D pupil detection confidence.
        """
        self.csv_export_filename = "blinks.csv"
        self.filter_len_s = filter_len_s
        self.onset_conf_thr = onset_conf_thr
        self.offset_conf_thr = offset_conf_thr

    @classmethod
    def dict_export(cls, raw_value: dict) -> dict:
        """
        获取字典值
        :param raw_value: 解包后的pupil datum数据
        :return: 待导出的字典值
        """
        blinks_id = raw_value["id"]
        timestamp = raw_value["timestamp"]
        confidence = raw_value["confidence"]
        return {"id": blinks_id, "timestamp": timestamp, "confidence": confidence}

    def process_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                           overwrite=False, n_jobs=None) -> typing.Union[pd.DataFrame, list]:
        """
        处理单个或多个瞳孔测量记录（多个时，可选并行或非并行）
        :param recordings: 处理单/多个记录，可以是包含瞳孔记录文件的单个文件夹，
                           也可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 处理后的所有数据 pd.DataFrame类型（单个瞳孔测量记录时）或list类型（多个瞳孔测量记录时）
        """
        if type(recordings) == list:
            res = self.process_multi_recordings(recordings, csv_out_dir, overwrite, n_jobs)
        else:
            rec_l = []
            for recording in os.listdir(recordings):
                if os.path.isfile(os.path.join(recordings, recording)):
                    rec_l = []
                    break
                rec_l.append(os.path.join(recordings, recording))
            if rec_l:
                res = self.process_multi_recordings(recordings, csv_out_dir, overwrite, n_jobs)
            else:
                res = self.process_single_recording(recordings, csv_out_dir, overwrite)
        return res

    def process_single_recording(self, recording: str, csv_out_dir: str = '', overwrite=False) -> typing.Optional[pd.DataFrame]:
        """
        处理单个瞳孔测量记录
        :param recording: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径或为带文件名的路径）
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :return: 所有数据 pd.DataFrame类型
        """
        if csv_out_dir == '':
            export_f = os.path.join(recording, 'exports', os.path.basename(recording), self.csv_export_filename)
        elif csv_out_dir.endswith('.csv'):
            export_f = csv_out_dir
        else:
            export_f = os.path.join(csv_out_dir, os.path.basename(recording.rstrip('/')), self.csv_export_filename)
        if os.path.exists(export_f):
            if not overwrite:
                print("{} exists already! Not overwriting, skipping.".format(export_f))
                return None
            else:
                logger.warning("{} exists already! Overwriting.".format(export_f))
        try:
            pupil_data = pd.DataFrame(self.load_and_yield_data(recording))
            if pupil_data.empty:
                raise ValueError("No pupil data found, pupil.pldata file may be empty, please check!")
            pupil_data.sort_values(by="timestamp", ignore_index=True, inplace=True)
            blink_data = self.blink_detection(recording, pupil_data)
            if not os.path.exists(os.path.dirname(export_f)):
                os.makedirs(os.path.dirname(export_f))
            blink_data.to_csv(export_f, index=False, encoding='utf-8-sig')
            blink_report = pd.DataFrame({'history_length_s': [self.filter_len_s],
                                         'onset_confidence_threshold': [self.onset_conf_thr],
                                         'offset_confidence_threshold': [self.offset_conf_thr],
                                         'blinks_exported': [len(blink_data)]})
            blink_report.to_csv(os.path.join(os.path.dirname(export_f), 'blink_report.csv'),
                                index=False, encoding='utf-8-sig')
            return blink_data
        except FileNotFoundError as e:
            print(f"文件缺失/格式有误，略过：{e}")
            return None

    def process_multi_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                                 overwrite=False, n_jobs=None) -> list:
        """
        处理多个瞳孔测量记录（并行或非并行）
        :param recordings: 处理多个记录，可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 并行处理后的所有数据 list类型
        """
        assert os.path.isdir(csv_out_dir) or csv_out_dir == '', "csv_out_dir仅接受已存在路径或空值''输入"
        assert (n_jobs is None) or (type(n_jobs) is int and n_jobs > 0) or (
                n_jobs == -1), 'n_jobs仅接受-1/正整数/None类型输入'
        if n_jobs == -1:
            n_jobs = None
        rec_l = []
        if type(recordings) == list:
            rec_l = recordings
        else:
            for recording in os.listdir(recordings):
                if os.path.isfile(os.path.join(recordings, recording)):
                    raise ValueError(u'输入文件夹内存在其他数据文件，请重新选择！')
                rec_l.append(os.path.join(recordings, recording))
        if n_jobs == 1:
            res = []
            for recording in rec_l:
                res.append(self.process_single_recording(recording, csv_out_dir, overwrite))
        else:
            parallel_process = partial(self.process_single_recording,
                                       csv_out_dir=csv_out_dir, overwrite=overwrite)
            with Pool(n_jobs) as pool:
                res = pool.map(parallel_process, rec_l)
        return res

    def load_and_yield_data(self, directory: str) -> typing.Generator:
        """
        Load and extract pupil diameter data
        :param directory: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :return: <class 'generator'> 所提取到的所有数据的生成器
        """
        msgpack_file = os.path.join(directory, "pupil.pldata")
        with open(msgpack_file, "rb") as fh:
            unpacker = msgpack.Unpacker(fh, raw=False, use_list=False)
            for topic, payload in unpacker:
                datum = msgpack.unpackb(payload, raw=False, use_list=False)
                if "2d" not in topic:
                    continue
                dict_row = self.dict_export(datum)
                yield dict_row

    def blink_detection(self, directory: str, pupil_data: pd.DataFrame) -> pd.DataFrame:
        """
        眨眼检测
        ref: https://github.com/pupil-labs/pupil/blob/master/pupil_src/shared_modules/blink_detection.py#L335-L458
        :param directory: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :param pupil_data: 瞳孔检测数据
        :return: 所有眨眼检测数据 pd.DataFrame类型
        """
        ts_start, ts_stop = pupil_data.timestamp.iloc[[0, -1]]
        total_time = ts_stop - ts_start
        activity = pupil_data.confidence
        filter_size = 2 * round(len(pupil_data) * self.filter_len_s / total_time / 2.0)
        blink_filter = np.ones(filter_size) / filter_size
        # This is different from the online filter. Convolution will flip the filter and
        # result in a reverse filter response. Therefore, we set the first half of the filter to -1 instead of
        # the second half such that we get the expected result.
        blink_filter[: filter_size // 2] *= -1
        # The theoretical response maximum is +-0.5
        # Response of +-0.45 seems sufficient for a confidence of 1.
        filter_response = fftconvolve(activity, blink_filter, "same") / 0.45
        onsets = filter_response > self.onset_conf_thr
        offsets = filter_response < -self.offset_conf_thr
        response_classification = np.zeros(filter_response.shape)
        response_classification[onsets] = 1.0
        response_classification[offsets] = -1.0
        # consolidation
        blink = None
        state = "no blink"  # others: 'blink started' | 'blink ending'
        blink_data = []
        counter = 1
        with open(os.path.join(directory, "info.player.json"), 'r') as f:
            info_data = json.load(f)
        try:
            world_ts = np.load(os.path.join(directory, "world_timestamps.npy"))
        except FileNotFoundError:
            logger.warning("world_timestamps.npy not exist already, "
                           "using the fake timestamps with 30Hz sampling rate of the world camera!")
            world_ts = np.arange(info_data['start_time_synced_s'],
                                 info_data['start_time_synced_s'] + info_data['duration_s'], 1/30.0)

        def blink_start(_idx):
            nonlocal blink, state, counter
            blink = {"start_timestamp": pupil_data.timestamp.iloc[_idx], "blink_id": counter, "__start_index__": _idx}
            state = "blink started"
            counter += 1

        def blink_finished(_idx):
            nonlocal blink
            start_index = blink["__start_index__"]
            blink["end_timestamp"] = pupil_data.timestamp.iloc[_idx]
            blink["duration_s"] = blink["end_timestamp"] - blink["start_timestamp"]
            blink["base_data"] = " ".join(map(str, pupil_data.timestamp[start_index:_idx]))
            blink["filter_response"] = " ".join(map(str, filter_response[start_index:_idx]))
            # blink confidence is the mean of the absolute filter response
            # during the blink event, clamped at 1.
            blink["confidence"] = min(float(np.abs(filter_response[start_index:_idx]).mean()), 1.0)
            # correlate world indices
            idx_start, idx_end = np.searchsorted(world_ts, [blink["start_timestamp"], blink["end_timestamp"]])
            idx_end = min(idx_end, len(world_ts) - 1)
            blink["start_frame_index"] = int(idx_start)
            blink["index"] = int((idx_start + idx_end) // 2)
            blink["end_frame_index"] = int(idx_end)
            blink_data.append(_BlinkDatum(id=blink["blink_id"],
                                          start_timestamp=blink["start_timestamp"],
                                          duration_s=blink["duration_s"],
                                          end_timestamp=blink["end_timestamp"],
                                          start_frame_index=blink["start_frame_index"],
                                          index=blink["index"],
                                          end_frame_index=blink["end_frame_index"],
                                          confidence=blink["confidence"],
                                          filter_response=blink["filter_response"],
                                          base_data=blink["base_data"]))
        for idx, classification in enumerate(response_classification):
            if state == "no blink" and classification > 0:
                blink_start(idx)
            elif state == "blink started" and classification == -1:
                state = "blink ending"
            elif state == "blink ending" and classification >= 0:
                blink_finished(idx - 1)  # blink ended previously
                if classification > 0:
                    blink_start(0)
                else:
                    blink = None
                    state = "no blink"
        if state == "blink ending":
            # only finish blink if it was already ending
            blink_finished(idx)  # idx is the last possible idx
        blink_data = pd.DataFrame(blink_data)
        return blink_data


class SurfaceDetectionExporter:
    """凝视和注视添加至表面检测数据提取
    由于这里的fixation等参数阈值未限制，默认为Pupil Player软件系统默认值，该项结果与Pupil Player提取的结果会有较大出入
    TODO: 待后续有空保证该结果一致"""
    def __init__(self):
        pass

    @classmethod
    def dict_export(cls, raw_value: dict, world_index: int, surface_topic: str) -> dict:
        """
        获取字典值
        :param raw_value: 解包后的pupil datum数据
        :param world_index: 世界视频帧索引
        :param surface_topic: surface主题名，gaze或fixations
        :return: 待导出的字典值
        """
        if surface_topic == "gaze":
            return {
                "world_timestamp": raw_value["world_timestamp"],
                "world_index": world_index,
                "gaze_timestamp": raw_value["timestamp"],
                "x_norm": raw_value["norm_pos"][0],
                "y_norm": raw_value["norm_pos"][1],
                "x_scaled": raw_value["x_scaled"],
                "y_scaled": raw_value["y_scaled"],
                "on_surf": raw_value["on_surf"],
                "confidence": raw_value["confidence"],
                "topic": raw_value["topic"],
            }
        else:  # fixations
            return {
                "world_timestamp": raw_value["world_timestamp"],
                "world_index": world_index,
                "fixation_id": raw_value["id"],
                "start_timestamp": raw_value["timestamp"],
                "duration": raw_value["duration"],
                "dispersion": raw_value["dispersion"],
                "norm_pos_x": raw_value["norm_pos"][0],
                "norm_pos_y": raw_value["norm_pos"][1],
                "x_scaled": raw_value["x_scaled"],
                "y_scaled": raw_value["y_scaled"],
                "on_surf": raw_value["on_surf"],
                "topic": raw_value["topic"],
            }

    def process_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                           overwrite=False, n_jobs=None) -> typing.Union[pd.DataFrame, list]:
        """
        处理单个或多个瞳孔测量记录（多个时，可选并行或非并行）
        :param recordings: 处理单/多个记录，可以是包含瞳孔记录文件的单个文件夹，
                           也可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 处理后的所有数据 pd.DataFrame类型（单个瞳孔测量记录时）或list类型（多个瞳孔测量记录时）
        """
        if type(recordings) == list:
            res = self.process_multi_recordings(recordings, csv_out_dir, overwrite, n_jobs)
        else:
            rec_l = []
            for recording in os.listdir(recordings):
                if os.path.isfile(os.path.join(recordings, recording)):
                    rec_l = []
                    break
                rec_l.append(os.path.join(recordings, recording))
            if rec_l:
                res = self.process_multi_recordings(recordings, csv_out_dir, overwrite, n_jobs)
            else:
                res = self.process_single_recording(recordings, csv_out_dir, overwrite)
        return res

    def process_single_recording(self, recording: str, csv_out_dir: str = '', overwrite=False) -> typing.Optional[dict]:
        """
        处理单个瞳孔测量记录
        :param recording: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径）
                            默认为对应recording路径/exports/recording文件名/surfaces文件夹
        :param overwrite: 是否对待保存的csv文件进行覆写
        :return: 所有数据 dict类型
        """
        if csv_out_dir == '':
            export_dir = os.path.join(recording, 'exports', os.path.basename(recording), 'surfaces')
        else:
            export_dir = os.path.join(csv_out_dir, os.path.basename(recording.rstrip('/')), 'surfaces')
        surface_data_dict = {}
        for surface_topic in "gaze", "fixations":
            try:
                surfaces_data = pd.DataFrame(self.load_and_yield_data(recording, surface_topic))
            except FileNotFoundError as e:
                print(f"文件缺失/格式有误，略过：{e}")
                return None
            if surface_topic == "gaze":
                surfaces_data.sort_values(by=["world_timestamp", "gaze_timestamp"], ignore_index=True, inplace=True)
            else:
                surfaces_data.sort_values(by=["world_timestamp", "start_timestamp"], ignore_index=True, inplace=True)
            if surfaces_data.empty:
                logger.warning(f"The recording {recording} contained no prerecorded "
                               f"{surface_topic}_on_surface, skiped.")
                continue
            for surface_name, surface_data in surfaces_data.groupby(surfaces_data.topic):
                surface_name = surface_name.split(".")[1]
                if not os.path.exists(export_dir):
                    os.makedirs(export_dir)
                export_f = os.path.join(export_dir, surface_topic + "_on_surface_" + surface_name + ".csv")
                if os.path.exists(export_f):
                    if not overwrite:
                        print("{} exists already! Not overwriting, skipping.".format(export_f))
                        return None
                    else:
                        logger.warning("{} exists already! Overwriting.".format(export_f))
                surface_data.drop("topic", axis=1, inplace=True)
                surface_data.to_csv(export_f, index=False, encoding='utf-8-sig')
                surface_data_dict[surface_name] = surface_data
        return surface_data_dict

    def process_multi_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                                 overwrite=False, n_jobs=None) -> list:
        """
        处理多个瞳孔测量记录（并行或非并行）
        :param recordings: 处理多个记录，可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 并行处理后的所有数据 list类型
        """
        assert os.path.isdir(csv_out_dir) or csv_out_dir == '', "csv_out_dir仅接受已存在路径或空值''输入"
        assert (n_jobs is None) or (type(n_jobs) is int and n_jobs > 0) or (
                n_jobs == -1), 'n_jobs仅接受-1/正整数/None类型输入'
        if n_jobs == -1:
            n_jobs = None
        rec_l = []
        if type(recordings) == list:
            rec_l = recordings
        else:
            for recording in os.listdir(recordings):
                if os.path.isfile(os.path.join(recordings, recording)):
                    raise ValueError(u'输入文件夹内存在其他数据文件，请重新选择！')
                rec_l.append(os.path.join(recordings, recording))
        if n_jobs == 1:
            res = []
            for recording in rec_l:
                res.append(self.process_single_recording(recording, csv_out_dir, overwrite))
        else:
            parallel_process = partial(self.process_single_recording,
                                       csv_out_dir=csv_out_dir, overwrite=overwrite)
            with Pool(n_jobs) as pool:
                res = pool.map(parallel_process, rec_l)
        return res

    def load_and_yield_data(self, directory: str, surface_topic: str) -> typing.Generator:
        """
        Load and extract pupil diameter data
        :param directory: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :param surface_topic: surface主题名，gaze或fixations
        :return: <class 'generator'> 所提取到的所有数据的生成器
        """
        data_ts = np.load(os.path.join(directory, "surfaces_timestamps.npy"))
        with open(os.path.join(directory, "info.player.json"), 'r') as f:
            info_data = json.load(f)
        try:
            world_ts = np.load(os.path.join(directory, "world_timestamps.npy"))
        except FileNotFoundError:
            logger.warning("world_timestamps.npy not exist already, "
                           "using the fake timestamps with 30Hz sampling rate of the world camera!")
            world_ts = np.arange(info_data['start_time_synced_s'],
                                 info_data['start_time_synced_s'] + info_data['duration_s'], 1/30.0)
        export_world_idc = find_closest(world_ts, data_ts)  # 获取最近的世界视频帧
        try:
            with open(os.path.join(directory, "surface_definitions_v01"), "rb") as fh:
                unpacker = msgpack.Unpacker(fh, raw=False, use_list=False)
                for payload in unpacker:
                    surf_def = payload["surfaces"]
        except FileNotFoundError:
            logger.warning("surface_definitions_v01 not exist already, using real_world_size = {'x': 1.0, 'y': 1.0}")
            surf_def = ({'name': '__DEFAULT__', 'real_world_size': {'x': 1.0, 'y': 1.0}}, )
        msgpack_file = os.path.join(directory, "surfaces.pldata")
        with open(msgpack_file, "rb") as fh:
            unpacker = msgpack.Unpacker(fh, raw=False, use_list=False)
            for ts, (topic, payload), idx in zip(data_ts, unpacker, export_world_idc):
                datum = msgpack.unpackb(payload, raw=False, use_list=False)
                surface_name = datum["name"]
                if surface_topic == "gaze":
                    for dat in datum["gaze_on_surfaces"]:
                        dat["world_timestamp"] = ts
                        dat["topic"] = topic
                        for surf in surf_def:
                            if surf["name"] in (surface_name, '__DEFAULT__'):
                                dat["x_scaled"] = dat["norm_pos"][0] * surf["real_world_size"]["x"]
                                dat["y_scaled"] = dat["norm_pos"][1] * surf["real_world_size"]["y"]
                        dict_row = self.dict_export(dat, idx, surface_topic)
                        yield dict_row
                else:  # fixations
                    for dat in datum["fixations_on_surfaces"]:
                        dat["world_timestamp"] = ts
                        dat["topic"] = topic
                        for surf in surf_def:
                            if surf["name"] in (surface_name, '__DEFAULT__'):
                                dat["x_scaled"] = dat["norm_pos"][0] * surf["real_world_size"]["x"]
                                dat["y_scaled"] = dat["norm_pos"][1] * surf["real_world_size"]["y"]
                        dict_row = self.dict_export(dat, idx, surface_topic)
                        yield dict_row


class PupilDataExporter:
    """Pupil数据提取"""
    def __init__(self, task: typing.Union[None, typing.List[str]] = None,
                 min_conf_thrs: typing.Union[None, typing.List[float]] = None,
                 filter_len_s: float = 0.20, onset_conf_thr: float = 0.50, offset_conf_thr: float = 0.50,
                 max_dispersion_deg: float = 1.50, min_duration_ms: int = 100, max_duration_ms: int = 400):
        """
        :param task: 待提取的任务，list类型，仅提取该列表中的任务对应的数据，默认为None，即提取全部任务的数据
        :param min_conf_thrs: 最小置信阈值列表，每个元素对应相应的任务的阈值，默认为None，取各自的默认值；
                      若列表长度小于任务数，则不够的任务阈值用最后一位阈值列表元素填充；若列表长度大于任务数，则截断
        :param filter_len_s: The time window's length in which the detector tries to find confidence drops and gains.
        :param onset_conf_thr: The threshold that the filter response ('Activity') must rise above to classify
                               the onset of a blink, corresponding to a sudden drop in 2D pupil detection confidence.
        :param offset_conf_thr: The threshold that the filter response ('Activity') must fall below to classify
                                the end of a blink, corresponding to a rise in 2D pupil detection confidence.
        :param max_dispersion_deg: Maximum distance between all gaze locations during a fixation.
        :param min_duration_ms: The minimum duration in which the dispersion threshold must not be exceeded.
        :param max_duration_ms: The maximum duration in which the dispersion threshold must not be exceeded.
        """
        tasks = ['pupil', 'gaze', 'blink', 'fixation', 'surface', ]
        if task is None or task == []:
            task = tasks
        assert set(task).issubset(set(tasks)), f'当前仅接受以下任务集中所列：{tasks}'
        if min_conf_thrs is None or min_conf_thrs == []:
            min_conf_thrs = [None] * len(task)
        elif len(min_conf_thrs) < len(task):
            min_conf_thrs = min_conf_thrs + [min_conf_thrs[-1]] * (len(task) - len(min_conf_thrs))
        else:
            min_conf_thrs = min_conf_thrs[:len(task)]
        min_conf_thr_dict = dict(zip(task, min_conf_thrs))
        self.data_exp = {}
        for i_t in task:
            if i_t == 'pupil':
                self.data_exp[i_t] = PupilDetectionExporter(min_conf_thr=min_conf_thr_dict[i_t])
            elif i_t == 'gaze':
                self.data_exp[i_t] = GazeDetectionExporter(min_conf_thr=min_conf_thr_dict[i_t])
            elif i_t == 'blink':
                blink = BlinkDetectionExporter(filter_len_s, onset_conf_thr, offset_conf_thr)
                self.data_exp[i_t] = blink
            elif i_t == 'fixation':
                fixation = FixationDetectionExporter(max_dispersion_deg, min_duration_ms,
                                                     max_duration_ms, min_conf_thr=min_conf_thr_dict[i_t])
                self.data_exp[i_t] = fixation
            elif i_t == 'surface':
                surface = SurfaceDetectionExporter()
                self.data_exp[i_t] = surface

    def process_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                           overwrite=False, n_jobs=None) -> dict:
        """
        处理单个或多个瞳孔测量记录（多个时，可选并行或非并行）
        :param recordings: 处理单/多个记录，可以是包含瞳孔记录文件的单个文件夹，
                           也可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 处理后的所有数据 dict类型
        """
        res_dict = {}
        for i_exp in self.data_exp:
            print(f'---------- Extracting {i_exp} data... ----------')
            res_dict[i_exp] = self.data_exp[i_exp].process_recordings(recordings, csv_out_dir, overwrite, n_jobs)
        return res_dict

    def process_single_recording(self, recording: str, csv_out_dir: str = '', overwrite=False) -> dict:
        """
        处理单个瞳孔测量记录
        :param recording: 处理单个记录，包含瞳孔记录文件的单个文件夹
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径或为带文件名的路径）
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :return: 所有数据 dict类型
        """
        res_dict = {}
        for i_exp in self.data_exp:
            print(f'---------- Extracting {i_exp} data... ----------')
            res_dict[i_exp] = self.data_exp[i_exp].process_single_recording(recording, csv_out_dir, overwrite)
        return res_dict

    def process_multi_recordings(self, recordings: typing.Union[str, list], csv_out_dir: str = '',
                                 overwrite=False, n_jobs=None) -> dict:
        """
        处理多个瞳孔测量记录（并行或非并行）
        :param recordings: 处理多个记录，可以是包含多个瞳孔记录文件的文件夹路径，也可以是包含多个瞳孔记录文件路径的列表
        :param csv_out_dir: 保存结果的CSV文件路径（为空或为路径），
                            默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件
        :param overwrite: 是否对待保存的csv文件进行覆写
        :param n_jobs: 并行运行CPU核数，默认为None;若为1非并行，若为-1或None,取os.cpu_count()全部核数,-1/正整数/None类型
        :return: 并行处理后的所有数据 dict类型
        """
        res_dict = {}
        for i_exp in self.data_exp:
            print(f'---------- Extracting {i_exp} data... ----------')
            res_dict[i_exp] = self.data_exp[i_exp].process_multi_recordings(recordings, csv_out_dir, overwrite, n_jobs)
        return res_dict


def find_closest(target, source):
    """Find indeces of closest `target` elements for elements in `source`.
    `source` is assumed to be sorted. Result has same shape as `source`.
    Implementation taken from:
    https://stackoverflow.com/questions/8914491/finding-the-nearest-value-and-return-the-index-of-array-in-python/8929827#8929827
    """
    target = np.asarray(target)  # fixes https://github.com/pupil-labs/pupil/issues/1439
    idx = np.searchsorted(target, source)
    idx = np.clip(idx, 1, len(target) - 1)
    left = target[idx - 1]
    right = target[idx]
    idx -= source - left < right - source
    return idx


def vector_dispersion(data):
    vectors = pd.DataFrame(data)[list("xyz")]
    distances = pdist(vectors, metric="cosine")
    dispersion = np.arccos(1.0 - distances.max())
    return np.rad2deg(dispersion)


def normalized2d_to_direction3d(point_2d, intrinsics):
    """
    Input shape: (2,)
    Output shape: (3,)
    """
    import cv2
    point_2d = np.asarray(point_2d, dtype="float64").reshape(2)
    camera_matrix = np.asarray(intrinsics["camera_matrix"])
    dist_coeffs = np.asarray(intrinsics["dist_coefs"])
    width, height = intrinsics["resolution"]
    # denormalize
    point_2d[:, 0] *= width
    point_2d[:, 1] = (1.0 - point_2d[:, 1]) * height
    # undistort
    point_2d = cv2.undistortPoints(point_2d, camera_matrix, dist_coeffs)
    # unproject
    point_3d = cv2.convertPointsToHomogeneous(point_2d)
    point_3d.shape = 3
    return point_3d


def direction3d_to_normalized2d(point_3d, intrinsics):
    """
    Input shape: (3,)
    Output shape: (2,)
    """
    import cv2
    point_3d = np.asarray(point_3d, dtype="float64").reshape((1, -1, 3))
    camera_matrix = np.asarray(intrinsics["camera_matrix"])
    dist_coeffs = np.asarray(intrinsics["dist_coefs"])
    width, height = intrinsics["resolution"]
    # rotation and translation of the camera, zero in our case
    rvec = tvec = np.zeros((1, 1, 3))
    # project and distort
    points_2d, _ = cv2.projectPoints(point_3d, rvec, tvec, camera_matrix, dist_coeffs)
    x, y = points_2d.reshape(2)
    # normalize
    x /= width
    y = 1.0 - y / height
    return x, y


def default_camera_intrinsics():
    # These are camera intrinsics that we recorded. They are estimates and generalize our
    # setup. Its always better to estimate intrinsics for each camera again.
    default_intrinsics = {
        "Pupil Cam1 ID2": {
            "version": 1,
            "(640, 480)": {
                "dist_coefs": [
                    [
                        -0.2430487205352619,
                        0.1623502095383119,
                        0.0001632500987373085,
                        8.322130878440475e-05,
                        0.017859803336754784,
                        0.1969284124154412,
                        0.00577741263771627,
                        0.09892258337410824,
                    ]
                ],
                "camera_matrix": [
                    [395.60662814306596, 0.0, 316.72212558212516],
                    [0.0, 395.56975615889445, 259.206579702132],
                    [0.0, 0.0, 1.0],
                ],
                'resolution': (640, 480),
                "cam_type": "radial",
            },
            "(1280, 720)": {
                "dist_coefs": [
                    [
                        -0.3758628065070806,
                        0.1643326166951343,
                        0.00012182540692089567,
                        0.00013422608638039466,
                        0.03343691733865076,
                        0.08235235770849726,
                        -0.08225804883227375,
                        0.14463365333602152,
                    ]
                ],
                "camera_matrix": [
                    [794.3311439869655, 0.0, 633.0104437728625],
                    [0.0, 793.5290139393004, 397.36927353414865],
                    [0.0, 0.0, 1.0],
                ],
                'resolution': (1280, 720),
                "cam_type": "radial",
            },
            "(1920, 1080)": {
                "dist_coefs": [
                    [
                        -0.13648546769272826,
                        -0.0033787366635030644,
                        -0.002343859061730869,
                        0.001926274947199097,
                    ]
                ],
                "camera_matrix": [
                    [793.8052697386686, 0.0, 953.2237035923064],
                    [0.0, 792.3104221704713, 572.5036513432223],
                    [0.0, 0.0, 1.0],
                ],
                'resolution': (1920, 1080),
                "cam_type": "fisheye",
            },
        },
        "Logitech Webcam C930e": {
            "(640, 480)": {
                "dist_coefs": [
                    [
                        0.10313391355051804,
                        -0.24657063652830105,
                        -0.001003806785350075,
                        -0.00046556297715377905,
                        0.1445780352338783,
                    ]
                ],
                "camera_matrix": [
                    [509.1810293948491, 0.0, 329.6996826114546],
                    [0.0, 489.7219438561515, 243.26037641451043],
                    [0.0, 0.0, 1.0],
                ],
                'resolution': (640, 480),
                "cam_type": "radial",
            },
            "(1280, 720)": {
                "dist_coefs": [
                    [
                        0.10152808562655541,
                        -0.23953332793667598,
                        -0.0021208895917640205,
                        -0.00023898995918166237,
                        0.1098748288957075,
                    ]
                ],
                "camera_matrix": [
                    [773.1676910077922, 0.0, 646.7114347564985],
                    [0.0, 743.1525324268981, 363.1646522363395],
                    [0.0, 0.0, 1.0],
                ],
                'resolution': (1280, 720),
                "cam_type": "radial",
            },
            "(1920, 1080)": {
                "dist_coefs": [
                    [
                        0.09961660299292627,
                        -0.21847900301383041,
                        -0.0010681464641609897,
                        -0.0014568525518904656,
                        0.09417837101183982,
                    ]
                ],
                "camera_matrix": [
                    [1120.4309938089518, 0.0, 968.3563459802797],
                    [0.0, 1077.3409390197398, 545.695766886239],
                    [0.0, 0.0, 1.0],
                ],
                'resolution': (1920, 1080),
                "cam_type": "radial",
            },
        },
        "PI world v1": {
            "(1088, 1080)": {
                "dist_coefs": [
                    [
                        -0.12390715699556255,
                        0.09983010007937897,
                        0.0013846287331131738,
                        -0.00036539454816030264,
                        0.020072404577046853,
                        0.2052173022520547,
                        0.009921380887245364,
                        0.06631870205961587,
                    ]
                ],
                "camera_matrix": [
                    [766.2927454396544, 0.0, 543.6272327745995],
                    [0.0, 766.3976103393867, 566.0580149497666],
                    [0.0, 0.0, 1.0],
                ],
                'resolution': (1088, 1080),
                "cam_type": "radial",
            }
        },
    }

    # Add measured intrinsics for the eyes (once for each ID for easy lookup)
    # TODO: From these intrinsics only the focal lengths were measured. The principal points
    # are just default values (half the resolution) and there's no distortion model for now.
    # At some later point we should measure the full intrinsics and replace existing
    # intrinsics with a recording upgrade.
    for eye_id in (0, 1):
        default_intrinsics.update(
            {
                f"Pupil Cam1 ID{eye_id}": {
                    "(320, 240)": {
                        "dist_coefs": [[0.0, 0.0, 0.0, 0.0, 0.0]],
                        "camera_matrix": [
                            [338.456035, 0.0, 160],
                            [0.0, 339.871543, 120],
                            [0.0, 0.0, 1.0],
                        ],
                        'resolution': (320, 240),
                        "cam_type": "radial",
                    },
                    "(640, 480)": {
                        "dist_coefs": [[0.0, 0.0, 0.0, 0.0, 0.0]],
                        "camera_matrix": [
                            [670.785555, 0.0, 320],
                            [0.0, 670.837798, 240],
                            [0.0, 0.0, 1.0],
                        ],
                        'resolution': (640, 480),
                        "cam_type": "radial",
                    },
                },
                f"Pupil Cam2 ID{eye_id}": {
                    "(192, 192)": {
                        "dist_coefs": [[0.0, 0.0, 0.0, 0.0, 0.0]],
                        "camera_matrix": [
                            [282.976877, 0.0, 96],
                            [0.0, 283.561467, 96],
                            [0.0, 0.0, 1.0],
                        ],
                        'resolution': (192, 192),
                        "cam_type": "radial",
                    },
                    "(400, 400)": {
                        "dist_coefs": [[0.0, 0.0, 0.0, 0.0, 0.0]],
                        "camera_matrix": [
                            [561.471804, 0.0, 200],
                            [0.0, 562.494105, 200],
                            [0.0, 0.0, 1.0],
                        ],
                        'resolution': (400, 400),
                        "cam_type": "radial",
                    },
                },
                f"Pupil Cam3 ID{eye_id}": {
                    "(192, 192)": {
                        "dist_coefs": [[0.0, 0.0, 0.0, 0.0, 0.0]],
                        "camera_matrix": [
                            [140.0, 0.0, 96],
                            [0.0, 140.0, 96],
                            [0.0, 0.0, 1.0],
                        ],
                        'resolution': (192, 192),
                        "cam_type": "radial",
                    },
                    "(400, 400)": {
                        "dist_coefs": [[0.0, 0.0, 0.0, 0.0, 0.0]],
                        "camera_matrix": [
                            [278.50, 0.0, 200],
                            [0.0, 278.55, 200],
                            [0.0, 0.0, 1.0],
                        ],
                        'resolution': (400, 400),
                        "cam_type": "radial",
                    },
                },
            }
        )
    return default_intrinsics


def parse_args():
    """
    命令行选项、参数和命令解析器
    :return: args_all，全部的参数解析
    """
    argparser = argparse.ArgumentParser(prog='PupilDataExporter',
                                        description=u'提取并保存Pupil数据. '
                                                    u'Same algorithm as the Pupil Player post-hoc blink detector.',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        usage='%(prog)s [-h] [-t TASKS] [-i RECORDINGS] '
                                              '[-o CSV_OUT_DIR] [-mct MIN_CONF_THRS] [-fl FILTER_LEN_S] '
                                              '[-onct ONSET_CONF_THR] [-offct OFFSET_CONF_THR] '
                                              '[-mdd MAX_DISPERSION_DEG] [-mind MIN_DURATION_MS] '
                                              '[-maxd MAX_DURATION_MS] [-w OVERWRITE] [-j N_JOBS]')
    argparser.add_argument('-t', '--tasks', metavar='', type=str, nargs='*', default=None,
                           help=u'待提取的任务（默认为None，即提取全部任务的数据）')
    argparser.add_argument('-i', '--recordings', metavar='', type=str,
                           help=u'输入单个或多个瞳孔记录文件')
    argparser.add_argument('-o', '--csv_out_dir', metavar='', type=str, default='',
                           help=u'保存结果的CSV文件路径'
                                u'（默认为对应recording路径/exports/recording文件名/self.csv_export_filename文件')
    argparser.add_argument('-mct', '--min_conf_thrs', metavar='', type=float, nargs='*', default=None,
                           help=u'最小置信阈值，每个元素对应相应的任务的阈值（默认为None，取各自的默认值）')
    argparser.add_argument('-fl', '--filter_len_s', metavar='', type=float, default=0.20,
                           help=u"The time window's length in which the detector tries to "
                                u"find confidence drops and gains.")
    argparser.add_argument('-onct', '--onset_conf_thr', metavar='', type=float, default=0.50,
                           help=u"The threshold that the filter response ('Activity') must rise above to classify the "
                                u"onset of a blink, corresponding to a sudden drop in 2D pupil detection confidence.")
    argparser.add_argument('-offct', '--offset_conf_thr', metavar='', type=float, default=0.50,
                           help=u"The threshold that the filter response ('Activity') must fall below to classify the "
                                u"end of a blink, corresponding to a rise in 2D pupil detection confidence.")
    argparser.add_argument('-mdd', '--max_dispersion_deg', metavar='', type=float, default=1.50,
                           help=u"Maximum distance between all gaze locations during a fixation.")
    argparser.add_argument('-mind', '--min_duration_ms', metavar='', type=int, default=80,
                           help=u"The minimum duration in which the dispersion threshold must not be exceeded.")
    argparser.add_argument('-maxd', '--max_duration_ms', metavar='', type=int, default=220,
                           help=u"The maximum duration in which the dispersion threshold must not be exceeded.")
    argparser.add_argument('-w', '--overwrite', metavar='', type=bool, default=False,
                           help=u'是否对结果进行覆写（默认为False）')
    argparser.add_argument('-j', '--n_jobs', metavar='', type=int, default=1,
                           help=u'并行核数（默认为1，不并行；-1为使用全部核心）')
    args_all = argparser.parse_args()
    argparser.print_help()  # 打印参数帮助信息
    return args_all


PUPIL_RAW_PATH = r'G:\Workspace\ZZK\tACS&MSIT&Stroop\PupilData'
PUPIL_EXPORTER_PATH = r'E:\pupil_anlys\anlys_res'

if __name__ == "__main__":
    # setup logging
    # logging.basicConfig(level=logging.DEBUG)
    # args = parse_args()
    # PupilDataExporter(args.tasks, args.min_conf_thrs, args.filter_len_s, args.onset_conf_thr,
    #                   args.offset_conf_thr, args.max_dispersion_deg, args.min_duration_ms,
    #                   args.max_duration_ms).process_recordings(args.recordings, args.csv_out_dir,
    #                                                            args.overwrite, args.n_jobs)
    # current_path = os.path.dirname(os.path.realpath(__file__))
    # pl_path = os.path.join(current_path, 'oculomotor_data')
    # PupilDataExporter([]).process_recordings(pl_path, csv_out_dir=r'./', overwrite=True, n_jobs=-1)

    def ext_raw(sub):
        pup_p = os.path.join(PUPIL_RAW_PATH, sub)
        exp_p = os.path.join(PUPIL_EXPORTER_PATH, sub)
        if not os.path.exists(exp_p):
            os.makedirs(exp_p)
        PupilDataExporter(['pupil']).process_recordings(pup_p, csv_out_dir=exp_p, overwrite=True, n_jobs=1)

    subj_l = [i for i in os.listdir(PUPIL_RAW_PATH) if os.path.isdir(os.path.join(PUPIL_RAW_PATH, i))]
    with Pool(os.cpu_count()) as pool:
        pool.map(ext_raw, subj_l)

