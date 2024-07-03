import argparse
import numpy as np
import random
import scipy
from . import utils
from . import make_param


def build_cylinder(args: argparse.Namespace):
    fr_list = args.resnames.split(":")
    fr_ratio = np.array(list(map(float,args.composition.split(":"))))/sum(list(map(float,args.composition.split(":"))))
    #fr周りの情報読み込み
    fr_atom_name_dict, fr_coord_dict, fr_resid_dict, fr_resname_dict = functional_residue.load_mol2(fr_list)
    #cylinderの生成
    cylinder_points,num_points,num_layer = cylinder.generate(args.radius, args.length, args.min_distance)
    #一層あたりのfunctional_residueの数
    fr_number_inlayer = functional_residue.get_number_inlayer(fr_ratio,num_points)
    output=""
    atom_index=1
    res_num=1
    for current_layer in range(1,num_layer+1):
        cylinder_points_layer = cylinder_points[num_points*(current_layer-1):num_points*current_layer]
        fr_order = functional_residue.get_order_inlayer(fr_list,fr_number_inlayer)
        #cylinderを1層ずつ修飾
        modified_fr_coord_list=cylinder.modify(cylinder_points_layer,fr_order,fr_coord_dict, outward=args.outward)
        for fr_name,modified_fr_coord in zip(fr_order,modified_fr_coord_list):
            atom_index=1
            for atom_name,atom_coord, resid, resname in zip(fr_atom_name_dict[fr_name],modified_fr_coord, fr_resid_dict[fr_name], fr_resname_dict[fr_name]):
                output+=(f"ATOM  {atom_index:5d}  {atom_name:4s}{resname}{resid:6d}    {atom_coord[0]:8.3f}{atom_coord[1]:8.3f}{atom_coord[2]:8.3f}  1.00  0.00\n")
                atom_index += 1
            output+="TER\n"
            res_num+=1
    np.savetxt(f"./cylinder.pdb",[output],fmt="%s")
    return None

def build_sphere(args: argparse.Namespace):
    fr_list = args.resnames.split(":")
    fr_ratio = np.array(list(map(float,args.composition.split(":"))))/sum(list(map(float,args.composition.split(":"))))
    #fr周りの情報読み込み
    fr_atom_name_dict, fr_coord_dict, fr_resid_dict, fr_resname_dict = functional_residue.load_mol2(fr_list)
    #cylinderの生成
    sphere_points,num_points = sphere.generate(args.radius,args.min_distance)
    output=""
    atom_index=1
    res_num=1
    fr_order=functional_residue.get_order_inlayer(fr_list,fr_ratio*100*num_points)
    sphere_center=np.mean(sphere_points, axis=0)
    for current_sphere_point,fr_name in zip(sphere_points,fr_order):
        fr_atom_name=fr_atom_name_dict[fr_name]
        fr_atom_coord=fr_coord_dict[fr_name]
        modify_coord=sphere.modify(current_sphere_point,fr_name,fr_atom_coord,sphere_center, outward=args.outward)
        for atom_name, atom_coord, resid, resname in zip(fr_atom_name, modify_coord, fr_resid_dict[fr_name], fr_resname_dict[fr_name]):
            output+=(f"ATOM  {atom_index:5d}  {atom_name:4s}{resname}{resid:6d}    {atom_coord[0]:8.3f}{atom_coord[1]:8.3f}{atom_coord[2]:8.3f}  1.00  0.00\n")
            atom_index+=1
        res_num+=1
        output+="TER\n"
    np.savetxt(f"./sphere.pdb",[output],fmt="%s")
    return None 

class calcurate_coordinate:
    def translation(modify_coord,fr_coord):
        C1_coord = fr_coord[0]
        translation_vector=np.array([modify_coord[0]-C1_coord[0],modify_coord[1]-C1_coord[1],modify_coord[2]-C1_coord[2]])
        translated_fr_coord = fr_coord + translation_vector
        return translated_fr_coord
    
    def rotation(coords: list, ref_axis: list): 
        original_fr_center= coords[0]
        coords = np.array(coords)-original_fr_center

        major_axis=coords[1]-coords[0]
        axis = np.cross(major_axis, ref_axis)
        axis_length = np.linalg.norm(axis)
        if axis_length == 0:
            return coords
        axis = axis/axis_length

        angle = np.arccos(np.dot(major_axis, ref_axis) / (np.linalg.norm(major_axis) * np.linalg.norm(ref_axis)))

        # Rodrigues' rotation formula
        K = np.array([[0, -axis[2], axis[1]],
                        [axis[2], 0, -axis[0]],
                        [-axis[1], axis[0], 0]])
        rotation_matrix = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

        rotated_fr_coord = np.dot(coords, rotation_matrix.T)

        # Random rotation to avoid atomic crash
        random_theta = np.random.uniform(0, np.pi)
        random_axis = ref_axis
        random_axis = random_axis / np.linalg.norm(random_axis)
        K = np.array([[0, -random_axis[2], random_axis[1]],
                        [random_axis[2], 0, -random_axis[0]],
                        [-random_axis[1], random_axis[0], 0]])
        rotation_matrix = np.eye(3) + np.sin(random_theta) * K + (1 - np.cos(random_theta)) * np.dot(K, K)
        rotated_fr_coord_randomized = np.dot(rotated_fr_coord, rotation_matrix.T)

        return rotated_fr_coord_randomized + original_fr_center
    
class cylinder:
    def generate(r, l, d):
        theta_spacing = d / r
        num_theta = int(np.ceil(2 * np.pi / theta_spacing))
        actual_theta_spacing = 2 * np.pi / num_theta

        z_spacing = np.sqrt(3) / 2 * d
        num_z = int(np.ceil(l / z_spacing))
        actual_z_spacing = l / num_z
        
        points = []
        for i in range(num_z):
            z = i * actual_z_spacing
            theta_offset = (i * actual_theta_spacing / 2) % (2 * np.pi)
            for j in range(num_theta):
                theta = (j * actual_theta_spacing + theta_offset) % (2 * np.pi)
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                points.append((x, y, z))
        return np.array(points)-np.mean(points, axis=0), num_theta, num_z
    
    def modify(cylinder_points_layer,fr_order,fr_coord_dict, outward=False):
        modify_coord_list=[]
        layer_center=np.average(cylinder_points_layer,axis=0)
        for modify_coord,fr_name in zip(cylinder_points_layer,fr_order):
            fr_coord = fr_coord_dict[fr_name]
            translated_fr_coord=calcurate_coordinate.translation(modify_coord,fr_coord)
            ref_vector = layer_center-modify_coord
            if outward == True:
                ref_vector = -1 * ref_vector
            rotated_fr_coord=calcurate_coordinate.rotation(translated_fr_coord,ref_vector)
            modify_coord_list.append(rotated_fr_coord)
        return modify_coord_list

class sphere:
    def estimate_number_of_points(r, d):
        sphere_area = 4 * np.pi * r ** 2
        area_per_point = d ** 2
        return int(np.ceil(sphere_area / area_per_point))

    def generate(r, d):
        num_points = sphere.estimate_number_of_points(r, d)
        points = []
        phi = np.pi * (3. - np.sqrt(5.))  # 黄金角

        for i in range(num_points):
            y = 1 - (i / float(num_points - 1)) * 2  # yは-1から1
            radius = np.sqrt(1 - y * y)  # yの高さにおける円の半径

            theta = phi * i  # 黄金角による角度

            x = np.cos(theta) * radius
            z = np.sin(theta) * radius

            points.append((r * x, r * y, r * z))
        
        return np.array(points)-np.mean(points, axis=0), num_points
    
    def modify(modify_coord,fr_name,fr_coord,sphere_center, outward=False):
        translated_fr_coord = calcurate_coordinate.translation(modify_coord,fr_coord)
        ref_vector = sphere_center-modify_coord
        if outward == True:
            ref_vector = -1 * ref_vector
        rotated_fr_coord = calcurate_coordinate.rotation(translated_fr_coord,ref_vector)
        modify_coord = rotated_fr_coord
        return modify_coord

class functional_residue:
    #1層に含まれるfunctional_residueの数を決める
    def get_number_inlayer(fr_ratio,num_points):
        fg_number = fr_ratio*100*num_points//num_points
        if sum(fg_number) != num_points:
            for i in np.argsort(np.remainder(fr_ratio*100*num_points,num_points))[::-1]:
                fg_number[i] += 1
                if sum(fg_number) == num_points:
                    break
        return np.array(fg_number)
    
    #１層に含まれるfunctional_residueの順番をrandamに決める
    def get_order_inlayer(fr_list,fg_number_inlayer):
            fr_order=[]
            for i,j in zip(fr_list,fg_number_inlayer):
                fr_order.extend([i]*int(j))
            random.shuffle(fr_order)
            return fr_order
    
    #functional_residueのmol2ファイルを読み込む
    def load_mol2_file(name):
        available_fr_list=make_param._get_params()
        
        if name not in available_fr_list:
            raise ValueError(f"{name} is not available")
        
        file_name=f"{available_fr_list[name]}/{name}.mol2"
        fr_coord = []
        fr_atom_name = []
        fr_resid = []
        fr_resname = []
        with open(file_name) as mol2_file:
            for atom_data in mol2_file:
                if atom_data.startswith('@<TRIPOS>ATOM'):
                    atom_section = True
                    continue
                elif atom_data.startswith('@<TRIPOS>'):
                    atom_section = False

                if atom_section:
                    atom_info = atom_data.split()
                    atom_name = atom_info[1]
                    coord = list(map(float,atom_info[2:5]))
                    fr_atom_name.append(atom_name)
                    fr_coord.append(coord)
                    fr_resname.append(atom_info[7])
                    fr_resid.append(int(atom_info[6]))
        ## 辞書にして簡素化できるはず
        return np.array(fr_atom_name),np.array(fr_coord), np.array(fr_resid), np.array(fr_resname)

    #functional_residueのmol2ファイルをlistとして受ける
    def load_mol2(fr_list):
        fr_atom_name_dict = {}
        fr_coord_dict = {}
        fr_resid_dict = {}
        fr_resname_dict = {}
        for fr_name in fr_list:
            fr_atom_name_dict[fr_name],fr_coord_dict[fr_name], fr_resid_dict[fr_name], fr_resname_dict[fr_name] = functional_residue.load_mol2_file(fr_name)
        return fr_atom_name_dict,fr_coord_dict, fr_resid_dict, fr_resname_dict