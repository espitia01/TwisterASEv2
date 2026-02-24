import numpy as np
import os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from scipy import spatial
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.path import Path
import matplotlib.patches as mpatches


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


def fold_into_cell(pos, A1, A2):
    """
    Fold Cartesian positions into the unit cell defined by A1 and A2.
    Returns positions with fractional parts in [0, 1) along A1 and A2.
    """
    cell2d = np.array([[A1[0], A2[0]],
                       [A1[1], A2[1]]])
    frac = np.linalg.solve(cell2d, pos[:, :2].T).T
    frac = frac % 1.0
    cart = frac[:, 0:1] * A1[:2] + frac[:, 1:2] * A2[:2]
    result = pos.copy()
    result[:, 0] = cart[:, 0]
    result[:, 1] = cart[:, 1]
    return result


def repeat_atoms(pos,rep,A1, A2, shift_origin):
    pos_sc = []
    origin = shift_origin[0]*A1 + shift_origin[1]*A2
    for i in range(-1,rep[0] - 1):
        for j in range(-1,rep[1] - 1):
            for p in pos:
                pos_sc.append(p + i*A1 + j*A2 - origin)
    return np.array(pos_sc)


def repeat_atoms_padded(pos, rep, A1, A2, shift_origin, padding=2):
    """
    Repeat atoms with padding in all directions to ensure plot is fully padded.
    
    Parameters:
    -----------
    pos : array
        Atom positions
    rep : tuple
        Number of repetitions in each direction
    A1, A2 : array
        Lattice vectors
    shift_origin : tuple
        Origin shift
    padding : int
        Number of extra unit cells to add in all directions for padding
    """
    pos_sc = []
    origin = shift_origin[0]*A1 + shift_origin[1]*A2
    
    for i in range(-padding, rep[0] + padding):
        for j in range(-padding, rep[1] + padding):
            for p in pos:
                pos_sc.append(p + i*A1 + j*A2 - origin)
    return np.array(pos_sc)


def draw_unit_cell(A1, A2, ax, color='black', linewidth=1.0, linestyle='-'):
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


def Rotate_one(a1, axis, theta):
    # Returns the rotation matrix associated with counterclockwise rotation about
    # a given axis by theta in radians.
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


def multiline(xs, ys, c, ax=None, **kwargs):
    # find axes
    ax = plt.gca() if ax is None else ax
    
    # create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)
    
    # set coloring of line segments
    lc.set_array(np.asarray(c))
    
    # add lines to axes and rescale
    ax.add_collection(lc)
    ax.autoscale()
    return lc


def strain_map(at_pos, ind, dist, alat):
    xs = []
    ys = []
    strain = []
    # Compute the local distance between each atom and the six nearest atoms
    for i in range(len(at_pos)):
        for j in range(len(ind[i])):
            # If the atom is not the same atom and the distance between the atom and the nearest atom is less than 4.0
            if ind[i,j] != i and dist[i,j] < 4.0:
                xs.append([at_pos[i,0], at_pos[ind[i,j],0]])
                ys.append([at_pos[i,1], at_pos[ind[i,j],1]])
                # Compute the strain in percentage change of alat
                strain.append((dist[i,j] - alat)*100/alat)
    return np.array(xs), np.array(ys), strain


def summarize_strain(strain, layer_name):
    s = np.array(strain)
    print(f"{layer_name} strain:   mean={s.mean():.2f}%,  med={np.median(s):.2f}%,  σ={s.std():.2f}%,  min={s.min():.2f}%,  max={s.max():.2f}%")
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    
    plt.figure()
    plt.hist(s, bins=50, edgecolor='k')
    plt.title(f"{layer_name} strain distribution", fontsize=14, family='Arial')
    plt.xlabel("Strain (%)", fontsize=12, family='Arial')
    plt.ylabel("Frequency", fontsize=12, family='Arial')
    plt.tight_layout()
    plt.savefig(f"hist_strain_{layer_name}.png", dpi=600, bbox_inches='tight')
    plt.savefig(f"hist_strain_{layer_name}.pdf", dpi=600, bbox_inches='tight')
    plt.close()


def _make_repeat_clip(A1, A2, repeat, ax):
    """
    Return a Polygon patch that covers exactly repeat[0]*A1 + repeat[1]*A2 cells.
    Used as a clip path so data outside the parallelogram is hidden.
    """
    O  = np.array([0.0, 0.0])
    P1 = repeat[0] * A1[:2]
    P2 = repeat[1] * A2[:2]
    P3 = P1 + P2
    verts = np.array([O, P1, P3, P2, O])
    return Polygon(verts, closed=True, transform=ax.transData)


def plot_strain_map(xs, ys, strain, A1, A2, xl, yl, place_dots, fnt_sz, tick_labsz, fig_name, repeat=None, show_unit_cell=False, show=False):
    """
    Publication-quality strain map plot with coolwarm colormap and optional unit cell.
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect('equal', 'box')
    
    lc = multiline(xs, ys, strain, cmap='coolwarm', lw=0.7)
    
    axcb = fig.colorbar(lc, ax=ax)
    axcb.set_label('Strain (%)', fontsize=fnt_sz, family='Arial')
    axcb.ax.tick_params(labelsize=tick_labsz, width=2.4, length=4.5)
    
    plt.xlim(xl[0], xl[1])
    plt.ylim(yl[0], yl[1])
    plt.xlabel(r'$x$ (Å)', fontsize=fnt_sz, family='Arial')
    plt.ylabel(r'$y$ (Å)', fontsize=fnt_sz, family='Arial')
    
    if show_unit_cell:
        draw_unit_cell(A1, A2, ax, color='black', linewidth=1.5)
    
    if place_dots:
        placedots(A1, A2)
    
    ax.tick_params(labelsize=tick_labsz, width=2.4, length=4.5)
    plt.tight_layout()
    plt.savefig(f"{fig_name}.png", format='png', dpi=600, bbox_inches='tight')
    plt.savefig(f"{fig_name}.pdf", format='pdf', dpi=600, bbox_inches='tight')
    plt.close()
    if show:
        plt.show()


def process_strain_directory(input_dir, show_unit_cell=False):
    """Process a single strain analysis directory."""
    print(f"\nProcessing strain analysis in: {input_dir}")
    
    original_dir = os.getcwd()
    os.chdir(input_dir)
    
    try:
        alat, A1, A2, A3, spc, rot_plot, rot_angle, xl, yl, place_dots, fnt_sz, tick_labsz, repeat, shift_origin = readinp_map('input')
        A1 = A1*alat
        A2 = A2*alat
        A3 = A3*alat
        
        pos = read_atoms('pos', spc)
        pos = np.array(pos)

        pos = fold_into_cell(pos, A1, A2)
        pos = repeat_atoms_padded(pos, repeat, A1, A2, shift_origin, padding=2)
        
        # Shift data so the visible region starts near origin
        # The visible region corner (0,0) in pre-shift coords is at -origin
        origin = shift_origin[0]*A1 + shift_origin[1]*A2
        shift_x = -origin[0]
        shift_y = -origin[1]
        pos[:, 0] -= shift_x
        pos[:, 1] -= shift_y
        
        # Compute visible region limits
        corners_x = []
        corners_y = []
        for ii in [0, repeat[0]]:
            for jj in [0, repeat[1]]:
                c = ii*A1 + jj*A2
                corners_x.append(c[0])
                corners_y.append(c[1])
        x_range = max(corners_x) - min(corners_x)
        y_range = max(corners_y) - min(corners_y)
        xl = [min(corners_x) - 0.05*x_range, max(corners_x) + 0.05*x_range]
        yl = [min(corners_y) - 0.05*y_range, max(corners_y) + 0.05*y_range]
        
        if rot_plot:
            axis = np.array([0,0,1])
            angle = rot_angle*np.pi/180.
            pos = Rotate_atoms(pos, axis, angle)
            A1 = Rotate_one(A1, axis, angle)
            A2 = Rotate_one(A2, axis, angle)
        
        pos_tree = spatial.cKDTree(pos)
        
        dist, ind = pos_tree.query(pos, k=7, distance_upper_bound=4)
        
        xs, ys, strain = strain_map(pos, ind, dist, alat)
        
        layer_name = os.path.basename(input_dir)
        
        summarize_strain(strain, layer_name)
        print(f"Max strain in {layer_name}: {np.max(strain):.2f}%")
        plot_strain_map(xs, ys, strain, A1, A2, xl, yl, place_dots, fnt_sz, tick_labsz, 
                       f"strain_{layer_name}", repeat=repeat, show_unit_cell=show_unit_cell, show=False)
        
        print(f"  - Generated strain plots in {input_dir}")
        
    except Exception as e:
        print(f"Error processing {input_dir}: {e}")
    finally:
        os.chdir(original_dir)


def main():
    """Main function to process all strain directories."""
    print("Processing strain analysis...")
    
    show_unit_cell = True
    
    base_path = '.'
    strain_path = os.path.join(base_path, 'StrainMap')
    
    if not os.path.exists(strain_path):
        print(f"Error: StrainMap directory not found at {strain_path}")
        print("Please run setup_analysis_dirs.py and generate_plot_inputs.py first.")
        return
    
    processed_count = 0
    for subdir in sorted(os.listdir(strain_path)):
        subdir_path = os.path.join(strain_path, subdir)
        if os.path.isdir(subdir_path):
            input_file = os.path.join(subdir_path, 'input')
            pos_file = os.path.join(subdir_path, 'pos')
            
            if all(os.path.exists(f) for f in [input_file, pos_file]):
                process_strain_directory(subdir_path, show_unit_cell)
                processed_count += 1
            else:
                print(f"Warning: Skipping {subdir_path} - missing required files")
    
    if processed_count == 0:
        print("No valid strain directories found to process.")
    else:
        print(f"\nCompleted processing {processed_count} strain directories.")


if __name__ == "__main__":
    main()
