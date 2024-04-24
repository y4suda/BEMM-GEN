import argparse
import numpy as np

from . import utils


def build_cylinder(args: argparse.Namespace):
    return None

def build_sphere(args: argparse.Namespace):
    return None 


def generate_cylinder(r, l, d):
    # 円周方向の角度間隔
    theta_spacing = d / r
    num_theta = int(np.ceil(2 * np.pi / theta_spacing))
    actual_theta_spacing = 2 * np.pi / num_theta
    
    # 層間の距離
    z_spacing = np.sqrt(3) / 2 * d
    num_z = int(np.ceil(l / z_spacing))
    actual_z_spacing = l / num_z
    
    points = []
    for i in range(num_z):
        z = i * actual_z_spacing
        # 各層で点を角度的にずらす（各層を1/2の角度ずらす）
        theta_offset = (i * actual_theta_spacing / 2) % (2 * np.pi)
        for j in range(num_theta):
            theta = (j * actual_theta_spacing + theta_offset) % (2 * np.pi)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            points.append((x, y, z))
    
    return np.array(points)