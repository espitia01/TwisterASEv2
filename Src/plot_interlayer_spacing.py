#!/usr/bin/python

import numpy as np
import os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from scipy import spatial
from mpl_toolkits.mplot3d import Axes3D


def readinp_map(file_name):
    fp = open(file_name,'r')
    lines = fp.readlines()
    fp.close()
    for i in range(len(lines)):
        # if line starts with #, ignore it
        if lines[i][0] == "#":
            continue
        if "alat" in lines[i]:
            w = lines[i].split("=") 
            alat = eval(w[1])
        if "A1" in lines[i]:
            w = lines[i].split("=")
            w = w[1].split()
            A1  = np.array([eval(w[0]), eval(w[1]), eval(w[2])])
        if "A2" in lines[i]:
            w = lines[i].split("=")
            w = w[1].split()
            A2  = np.array([eval(w[0]), eval(w[1]), eval(w[2])])
        if "A3" in lines[i]:
            w = lines[i].split("=")
            w = w[1].split()
            A3  = np.array([eval(w[0]), eval(w[1]), eval(w[2])])
        if "compute_spc" in lines[i]:
            w = lines[i].split("=")
            # Handle both [W W] and W W formats
            spc_str = w[1].strip()
            if '[' in spc_str and ']' in spc_str:
                # Remove brackets and split
                spc_str = spc_str.replace('[', '').replace(']', '')
            spc = spc_str.split()
        if "rotate_plot" in lines[i]:
            w = lines[i].split("=")
            rot_plot = eval(w[1])
        if "rotate_angle" in lines[i]:
            w = lines[i].split("=")
            rot_angle = eval(w[1])
        if "xlim" in lines[i]:
            w = lines[i].split("=")
            w = w[1].split()
            xl = [eval(w[0]),eval(w[1])]
        if "ylim" in lines[i]:
            w = lines[i].split("=")
            w = w[1].split()
            yl = [eval(w[0]), eval(w[1])]
        if "place_dots" in lines[i]:
            w = lines[i].split("=")
            place_dots = eval(w[1])
        if "fontsize" in lines[i]:
            w = lines[i].split("=")
            fnt_sz = eval(w[1])
        if "tick_labelsize" in lines[i]:
            w = lines[i].split("=")
            tick_labsz = eval(w[1])
        if "repeat_units" in lines[i]:
            w = lines[i].split("=")
            w = w[1].split()
            repeat = [eval(w[0]), eval(w[1])]
        if "shift_origin" in lines[i]:
            w = lines[i].split("=")
            w = w[1].split()
            shift_origin = [eval(w[0]), eval(w[1])]
    return alat, A1, A2, A3, spc, rot_plot, rot_angle, xl, yl, place_dots, fnt_sz, tick_labsz, repeat, shift_origin


def read_atoms(file_name, spc):
    fp = open(file_name,'r')
    lines = fp.readlines()
    fp.close()
    A = []
    for i in range(len(lines)):
        w = lines[i].split()
        if len(w) >= 4 and w[0] in spc:
            A.append([eval(w[1]),eval(w[2]),eval(w[3])])
    if len(A) == 0:
        print(f"Warning: No atoms of species {spc} found in {file_name}")
    return np.array(A)


def placedots(A1, A2, size = 60):
    plt.scatter(0 , 0 , s = size, color = "dodgerblue")
    ot = (A1 + A2)/3.
    plt.scatter(ot[0] , ot[1] , s = size, color = "red")
    tt = (A1 + A2)*2/3.
    plt.scatter(tt[0] , tt[1] , s = size, color = "green")
    Sum = A1 + A2
    plt.scatter(Sum[0], Sum[1], s = size, color = "dodgerblue")


def placeaxes(A1, A2):
    plt.plot([0 ,A1[0] ], [0 , A1[1]], color = "grey", lw = 2.1)
    plt.plot([0, A2[0]], [0,  A2[1]], color = "grey", lw = 2.1)
    Sum = A1 + A2
    plt.plot([ A2[0], Sum[0]  ], [A2[1], Sum[1]], color = "grey", lw = 2.1)
    plt.plot([ A1[0], Sum[0]], [A1[1], Sum[1]], color = "grey", lw = 2.1)


def repeat_data(w_xy,rep,A1,A2, shift_origin):
    w_xy_sc = []
    origin = -shift_origin[0]*A1 - shift_origin[1]*A2
    for i in range(-1,rep[0]-1):
        for j in range(-1,rep[1]-1):
            for k in range(len(w_xy)):
                xyr = np.array([w_xy[k][0],w_xy[k][1], 30.48780487804878])
                xyr = xyr + i*A1 + j*A2 - origin
                w_xy_sc.append([xyr[0],xyr[1],w_xy[k][2]])
    return np.array(w_xy_sc)


def repeat_data_padded(w_xy, rep, A1, A2, shift_origin, padding=2):
    """
    Repeat data with padding in all directions to ensure plot is fully padded.
    
    Parameters:
    -----------
    w_xy : array
        Data points with x, y, and color values
    rep : tuple
        Number of repetitions in each direction
    A1, A2 : array
        Lattice vectors
    shift_origin : tuple
        Origin shift
    padding : int
        Number of extra unit cells to add in all directions for padding
    """
    w_xy_sc = []
    origin = shift_origin[0]*A1 + shift_origin[1]*A2
    
    for i in range(-padding, rep[0] + padding):
        for j in range(-padding, rep[1] + padding):
            for k in range(len(w_xy)):
                xyr = np.array([w_xy[k][0], w_xy[k][1], 0.0])
                xyr = xyr + i*A1 + j*A2 - origin
                w_xy_sc.append([xyr[0], xyr[1], w_xy[k][2]])
    return np.array(w_xy_sc)


def draw_unit_cell(A1, A2, ax, color='black', linewidth=1.5, linestyle='-'):
    """
    Draw the unit cell boundaries.
    
    Parameters:
    -----------
    A1, A2 : array
        Lattice vectors
    ax : matplotlib axis
        Axis to draw on
    color : str
        Line color
    linewidth : float
        Line width
    linestyle : str
        Line style
    """
    ax.plot([0, A1[0]], [0, A1[1]], color=color, lw=linewidth, linestyle=linestyle, zorder=10)
    ax.plot([0, A2[0]], [0, A2[1]], color=color, lw=linewidth, linestyle=linestyle, zorder=10)
    Sum = A1 + A2
    ax.plot([A2[0], Sum[0]], [A2[1], Sum[1]], color=color, lw=linewidth, linestyle=linestyle, zorder=10)
    ax.plot([A1[0], Sum[0]], [A1[1], Sum[1]], color=color, lw=linewidth, linestyle=linestyle, zorder=10)


def create_blue_red_white_colormap():
    """
    Create a publication-quality blue-red-white diverging colormap.
    Blue (low) -> White (middle) -> Red (high)
    """
    colors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7', 
              '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list('blue_red_white', colors, N=n_bins)
    return cmap


def repeat_atoms(pos,rep,A1, A2, shift_origin):
    pos_sc = []
    origin = shift_origin[0]*A1 + shift_origin[1]*A2
    for i in range(-1,rep[0] - 1):
        for j in range(-1,rep[1] - 1):
            for p in pos:
                pos_sc.append(p + i*A1 + j*A2 - origin)
    return np.array(pos_sc) 


def process_interlayer_directory(input_dir, show_unit_cell=False):
    """Process a single interlayer spacing directory."""
    print(f"\nProcessing interlayer spacing analysis in: {input_dir}")
    
    original_dir = os.getcwd()
    os.chdir(input_dir)
    
    try:
        alat, A1, A2, A3, spc, rot_plot, rot_angle, xl, yl, place_dots, fnt_sz, tick_labsz, repeat, shift_origin = readinp_map('input')
        A1 = A1*alat
        A2 = A2*alat
        A3 = A3*alat
        
        pos_bot = read_atoms('pos_l', spc)
        pos_top = read_atoms('pos_u', spc)
        
        if rot_plot:
            axis = np.array([0,0,1])
            angle = rot_angle*np.pi/180.
            pos_bot = Rotate_atoms(pos_bot, axis, angle)
            pos_top = Rotate_atoms(pos_top, axis, angle)
            A1 = Rotate_one(A1, axis, angle)
            A2 = Rotate_one(A2, axis, angle)
        
        pos_top_tree = spatial.cKDTree(pos_top)
        dist, ind = pos_top_tree.query(pos_bot)
        
        d = np.zeros(len(pos_bot))
        for i in range(len(pos_bot)):
            d[i] = pos_top[ind[i],2] - pos_bot[i,2]
        print("min, max, mean:", np.amin(d), np.amax(d), np.mean(d))
        
        arr_xyc = np.zeros((len(pos_bot),3))
        arr_xyc[:,0] = pos_bot[:,0]
        arr_xyc[:,1] = pos_bot[:,1]
        arr_xyc[:,2] = d
        
        arr_xy_sc = repeat_data_padded(arr_xyc, repeat, A1, A2, shift_origin, padding=2)
        
        # Shift data so minimum position is at origin
        min_x = np.min(arr_xy_sc[:, 0])
        min_y = np.min(arr_xy_sc[:, 1])
        arr_xy_sc[:, 0] -= min_x
        arr_xy_sc[:, 1] -= min_y
        
        # Adjust plot limits accordingly
        xl = [xl[0] - min_x, xl[1] - min_x]
        yl = [yl[0] - min_y, yl[1] - min_y]
        
        layer_name = os.path.basename(input_dir)
        plot_scatter_publication(arr_xy_sc, A1, A2, xl, yl, place_dots, fnt_sz, tick_labsz, 
                                show_unit_cell, layer_name, show=False)
        plot_interpolated_publication(arr_xy_sc, A1, A2, xl, yl, place_dots, fnt_sz, tick_labsz, 
                                     show_unit_cell, layer_name, show=False)
        
        print(f"  - Generated plots in {input_dir}")
        
    except Exception as e:
        print(f"Error processing {input_dir}: {e}")
    finally:
        os.chdir(original_dir)


def main():
    """Main function to process all interlayer spacing directories."""
    print("Processing interlayer spacing analysis...")
    
    show_unit_cell = False
    
    base_path = '.'
    ils_path = os.path.join(base_path, 'InterlayerSpacingMap')
    
    if not os.path.exists(ils_path):
        print(f"Error: InterlayerSpacingMap directory not found at {ils_path}")
        print("Please run setup_analysis_dirs.py and generate_plot_inputs.py first.")
        return
    
    processed_count = 0
    for subdir in sorted(os.listdir(ils_path)):
        subdir_path = os.path.join(ils_path, subdir)
        if os.path.isdir(subdir_path):
            input_file = os.path.join(subdir_path, 'input')
            pos_l_file = os.path.join(subdir_path, 'pos_l')
            pos_u_file = os.path.join(subdir_path, 'pos_u')
            
            if all(os.path.exists(f) for f in [input_file, pos_l_file, pos_u_file]):
                process_interlayer_directory(subdir_path, show_unit_cell)
                processed_count += 1
            else:
                print(f"Warning: Skipping {subdir_path} - missing required files")
    
    if processed_count == 0:
        print("No valid interlayer spacing directories found to process.")
    else:
        print(f"\nCompleted processing {processed_count} interlayer spacing directories.")


def interp(q_xyc, xmin, xmax, ymin, ymax,dxi):
    X = q_xyc[:,0]
    Y = q_xyc[:,1]
    C = q_xyc[:,2]
    xi = np.arange(xmin,xmax,dxi)
    yi = np.arange(ymin,ymax,dxi)
    xi,yi = np.meshgrid(xi,yi)
    zi = griddata((X,Y),C,(xi,yi),method='linear',fill_value = 0)
    return xi, yi, zi


def Rotate_one( a1, axis, theta):
  #   Returns the rotation matrix associated with counterclockwise rotation about
  #   a given axis by  theta in radians.
  axis = np.asarray(axis)
  theta = np.asarray(theta)
  axis = axis/np.sqrt(np.dot(axis, axis))
  a = np.cos(theta/2)
  b, c, d = -axis*np.sin(theta/2)
  aa, bb, cc, dd = a*a, b*b, c*c, d*d
  bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
  M = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                   [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                   [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
  return np.dot(M,a1)


def Rotate_atoms(layer,norm,angle):
  layer_r = []
  for i in range(len(layer)):
    layer_r.append(Rotate_one( layer[i], norm, angle))
  return np.array(layer_r)


def d_fromline(p1,p2,p3):
  d = np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1)
  return d


def rep_cell(pos,a_1, a_2):
    a_1 = np.array([a_1[0],a_1[1],0.0])
    a_2 = np.array([a_2[0],a_2[1],0.0])
    pos_sc = []
    for i in range(-1,2):
        for j in range(-1,2):
            for p in pos:
                pos_sc.append(p + i*a_1 + j*a_2)  
    return np.array(pos_sc)



def plot_interpolated_publication(arr_xyc, A1, A2, xl, yl, place_dots, fnt_sz, tick_labsz, show_unit_cell, layer_name, show=False):
    """
    Publication-quality interpolated plot with coolwarm colormap and unit cell.
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    
    xmin, xmax = xl[0], xl[1]
    ymin, ymax = yl[0], yl[1]
    vmin, vmax = np.amin(arr_xyc[:, 2]), np.amax(arr_xyc[:, 2])
    
    xint, yint, zint = interp(arr_xyc, xmin, xmax, ymin, ymax, 0.15)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect('equal', 'box')
    
    hm = plt.pcolormesh(xint, yint, zint, cmap=cm.coolwarm, rasterized=True, 
                        shading='auto', vmin=vmin, vmax=vmax)
    
    if show_unit_cell:
        draw_unit_cell(A1, A2, ax, color='black', linewidth=1.5)
    
    if place_dots:
        placedots(A1, A2)
    
    cbar = plt.colorbar(hm, ax=ax, shrink=1.0, pad=0.02, aspect=20)
    cbar.set_label('Interlayer spacing (Å)', fontsize=fnt_sz, family='Arial')
    cbar.ax.tick_params(labelsize=tick_labsz, width=2.4, length=4.5)
    
    plt.xlabel(r'$x$ (Å)', fontsize=fnt_sz, family='Arial')
    plt.ylabel(r'$y$ (Å)', fontsize=fnt_sz, family='Arial')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    ax.tick_params(labelsize=tick_labsz, width=2.4, length=4.5)
    
    plt.tight_layout()
    plt.savefig(f'interpolated_{layer_name}.png', format='png', dpi=600, bbox_inches='tight')
    plt.savefig(f'interpolated_{layer_name}.pdf', format='pdf', dpi=600, bbox_inches='tight')
    plt.close()
    if show:
        plt.show()


def plot_interpolated(arr_xyc, A1, A2, xl, yl, place_dots, fnt_sz, tick_labsz, show = False):
    xmin, xmax = xl[0], xl[1]
    ymin, ymax = yl[0], yl[1]
    vmin, vmax = np.amin(arr_xyc[:,2]), np.amax(arr_xyc[:,2])
    xint, yint, zint = interp(arr_xyc, xmin, xmax, ymin, ymax,0.15)
    fig, ax = plt.subplots()
    hm = plt.pcolormesh(xint,yint,zint ,cmap = cm.inferno, rasterized=True, shading = 'auto', vmin = vmin, vmax = vmax)
    placeaxes(A1, A2)
    if place_dots:
        placedots(A1, A2)
    cbar = plt.colorbar(hm)
    cbar.ax.tick_params(labelsize=tick_labsz,width = 2.4, length=4.5)
    plt.xlabel(r"$x (\mathrm{\AA})$", fontsize = fnt_sz)
    plt.ylabel(r"$y (\mathrm{\AA}$)", fontsize = fnt_sz)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax.tick_params(labelsize = tick_labsz)
    plt.tight_layout()
    plt.savefig('interpolated.png',format = 'png', dpi = 400)
    if show:
        plt.show()




def plot_scatter_publication(arr_xy_sc, A1, A2, xl, yl, place_dots, fnt_sz, tick_labsz, show_unit_cell, layer_name, show=False):
    """
    Publication-quality scatter plot with coolwarm colormap and unit cell.
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect('equal', 'box')
    
    scatter = plt.scatter(arr_xy_sc[:, 0], arr_xy_sc[:, 1], s=5, 
                         c=arr_xy_sc[:, 2], cmap=cm.coolwarm, edgecolor='None')
    
    plt.xlim(xl[0], xl[1])
    plt.ylim(yl[0], yl[1])
    
    cbar = plt.colorbar(scatter, ax=ax, shrink=1.0, pad=0.02, aspect=20)
    cbar.set_label('Interlayer spacing (Å)', fontsize=fnt_sz, family='Arial')
    cbar.ax.tick_params(labelsize=tick_labsz, width=2.4, length=4.5)
    
    if show_unit_cell:
        draw_unit_cell(A1, A2, ax, color='black', linewidth=1.5)
    
    if place_dots:
        placedots(A1, A2, size=60)
    
    dot_xy = np.array([
        [0.0,          0.0],
        [A1[0],        A1[1]],
        [A2[0],        A2[1]],
        [A1[0]+A2[0],  A1[1]+A2[1]],
    ])
    d_interp = griddata(
        (arr_xy_sc[:,0], arr_xy_sc[:,1]),
        arr_xy_sc[:,2],
        dot_xy,
        method='cubic'
    )
    print("Interlayer distances at lattice-dot positions (interpolated):")
    for (xd, yd), dval in zip(dot_xy, d_interp):
        print(f"  • dot at ({xd:.2f}, {yd:.2f}) Å → d = {dval:.3f} Å")
    
    plt.xlabel(r'$x$ (Å)', fontsize=fnt_sz, family='Arial')
    plt.ylabel(r'$y$ (Å)', fontsize=fnt_sz, family='Arial')
    ax.tick_params(labelsize=tick_labsz, width=2.4, length=4.5)
    
    plt.tight_layout()
    plt.savefig(f'scatter_{layer_name}.png', format='png', dpi=600, bbox_inches='tight')
    plt.savefig(f'scatter_{layer_name}.pdf', format='pdf', dpi=600, bbox_inches='tight')
    plt.close()
    if show:
        plt.show()


def plot_scatter(arr_xy_sc, A1, A2, xl, yl, place_dots, fnt_sz, tick_labsz, show=False):
    fig, ax = plt.subplots()

    sc = ax.scatter(
        np.abs(arr_xy_sc[:,0]),
        np.abs(arr_xy_sc[:,1]),
        s=17,
        c=np.abs(arr_xy_sc[:,2]),
        edgecolor="None",
        cmap=cm.coolwarm
    )
    plt.xlim(xl[0], xl[1])
    plt.ylim(yl[0], yl[1])
    cbar = plt.colorbar(sc)
    cbar.ax.tick_params(labelsize=tick_labsz, width=2., length=4.5)
    ax.tick_params(labelsize=tick_labsz)

    placeaxes(A1, A2)
    if place_dots:
        placedots(A1, A2, size=60)

    dot_xy = np.array([
        [0.0,          0.0],
        [A1[0],        A1[1]],
        [A2[0],        A2[1]],
        [A1[0]+A2[0],  A1[1]+A2[1]],
    ])
    d_interp = griddata(
        (arr_xy_sc[:,0], arr_xy_sc[:,1]),
        arr_xy_sc[:,2],
        dot_xy,
        method='cubic'
    )
    print("Interlayer distances at your lattice‑dot positions (interpolated):")
    for (xd, yd), dval in zip(dot_xy, d_interp):
        print(f"  • dot at ({xd:.2f}, {yd:.2f}) Å → d = {dval:.3f} Å")

    plt.xlabel(r"$x\;(\mathrm{\AA})$", fontsize=fnt_sz)
    plt.ylabel(r"$y\;(\mathrm{\AA})$", fontsize=fnt_sz)
    plt.tight_layout()
    plt.savefig('scatter.png', format='png', dpi=400)
    if show:
        plt.show()





if __name__ == "__main__":
    main()
