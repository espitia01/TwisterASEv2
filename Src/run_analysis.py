"""
run_analysis.py - Complete analysis workflow for layered materials

Orchestrates the full analysis pipeline:
1. Setup analysis directories (InterlayerSpacingMap, StrainMap)
2. Generate input files for each analysis
3. Run interlayer spacing analysis (for multi-layer systems)
4. Run strain analysis (for all layers)
"""

import subprocess
import sys
import os
import time
import glob
from datetime import datetime


def run_script(script_name, description):
    """Run a Python script and report results."""
    print(f"\n{'='*60}")
    print(f"STEP: {description}")
    print(f"Running: {script_name}")
    print(f"{'='*60}")
    
    start_time = time.time()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(script_dir, script_name)
    
    try:
        result = subprocess.run([sys.executable, script_path], 
                              capture_output=True, 
                              text=True, 
                              check=True)
        
        if result.stdout:
            print(result.stdout)
        
        elapsed = time.time() - start_time
        print(f"[OK] {description} completed successfully in {elapsed:.2f}s")
        return True
        
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"[FAIL] {description} failed after {elapsed:.2f}s")
        print(f"Error code: {e.returncode}")
        if e.stdout:
            print("STDOUT:")
            print(e.stdout)
        if e.stderr:
            print("STDERR:")
            print(e.stderr)
        return False
    except FileNotFoundError:
        print(f"[FAIL] Script not found: {script_name}")
        return False


def check_prerequisites():
    """Check for required input files."""
    # Check for cutpos.inp
    if not os.path.exists('cutpos.inp'):
        print("[FAIL] Missing required file: cutpos.inp")
        return False
    
    # Check for layer input files
    layer_inp_files = glob.glob('layer*.inp')
    if len(layer_inp_files) < 1:
        print(f"[FAIL] Need at least 1 layer*.inp file, found {len(layer_inp_files)}")
        return False
    
    # Check for position files
    pos_files = glob.glob('layer_*_coords.dat')
    
    if len(pos_files) < 1:
        print(f"[FAIL] Need at least 1 layer_*_coords.dat file, found {len(pos_files)}")
        print("Run cutpos.py first to extract layer coordinates")
        return False
    
    print(f"[OK] Found prerequisite files:")
    print(f"   - cutpos.inp")
    print(f"   - {len(layer_inp_files)} layer input files: {layer_inp_files}")
    print(f"   - {len(pos_files)} position files: {pos_files}")
    
    # Check system type
    if len(pos_files) == 1:
        print("   Detected: Monolayer system (strain analysis only)")
    else:
        print(f"   Detected: {len(pos_files)}-layer system (strain + interlayer spacing analysis)")
    
    return True


def main():
    """Main function to run the complete analysis workflow."""
    print("LAYERED MATERIAL ANALYSIS WORKFLOW")
    print("=" * 60)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Working directory: {os.getcwd()}")
    
    # Check prerequisites
    if not check_prerequisites():
        sys.exit(1)
    
    # Define the analysis workflow
    pos_files = glob.glob('layer_*_coords.dat')
    
    workflow_steps = [
        ("setup_analysis_dirs.py", "Setting up analysis directories"),
        ("generate_plot_inputs.py", "Generating input files for analysis")
    ]
    
    # Only add interlayer spacing analysis for multi-layer systems
    if len(pos_files) > 1:
        workflow_steps.append(("plot_interlayer_spacing.py", "Running interlayer spacing analysis"))
    
    # Always add strain analysis
    workflow_steps.append(("plot_strains.py", "Running strain analysis"))
    
    # Track success/failure
    successful_steps = 0
    total_steps = len(workflow_steps)
    start_time = time.time()
    
    # Execute each step
    for script, description in workflow_steps:
        success = run_script(script, description)
        if success:
            successful_steps += 1
        else:
            print(f"\n[FAIL] Workflow stopped due to failure in: {description}")
            break
    
    # Summary
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print("WORKFLOW SUMMARY")
    print(f"{'='*60}")
    print(f"Completed: {successful_steps}/{total_steps} steps")
    print(f"Total time: {total_time:.2f}s")
    print(f"Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if successful_steps == total_steps:
        print("All analysis steps completed successfully!")
        print("\nGenerated outputs:")
        
        # List interlayer spacing outputs
        ils_dirs = glob.glob("InterlayerSpacingMap/Layer_*")
        if ils_dirs:
            print("  InterlayerSpacingMap/")
            for ils_dir in sorted(ils_dirs):
                layer_name = os.path.basename(ils_dir)
                print(f"     {layer_name}/")
                print(f"       - scatter.png (scatter plot)")
                print(f"       - interpolated.png (heatmap)")
        
        # List strain outputs
        strain_dirs = glob.glob("StrainMap/Layer_*")
        if strain_dirs:
            print("  StrainMap/")
            for strain_dir in sorted(strain_dirs):
                layer_name = os.path.basename(strain_dir)
                print(f"     {layer_name}/")
                print(f"       - strain_{layer_name}.png (strain map)")
                print(f"       - hist_strain_{layer_name}.png (strain histogram)")
        return 0
    else:
        print("[FAIL] Workflow incomplete due to errors")
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
