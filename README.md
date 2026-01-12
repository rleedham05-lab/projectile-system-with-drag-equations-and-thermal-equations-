# projectile-system-with-drag-equations-and-thermal-equations-
# importing packages 
import math 
import time 
import os 
import sys 
import numpy as np
import tkinter as tk
from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import messagebox 
from math import sqrt  
import matplotlib.pyplot as plt

# Global variables for entry widgets
entry_velocity = None
entry_mass = None
entry_angle = None
entry_emistivity = None
entry_specific_heat_capacity = None
entry_area = None
entry_density = None
entry_drag_coeffecient = None
entry_enviroment_temperature = None

# Global variables for labels
velocity_Label = None
mass_Label = None
angle_Label = None
emistivity_Label = None
specific_heat_capacity_Label = None
aera_Label = None
density_Label = None
drag_coeffecient = None
enviroment_temperature = None
Fdrag_Label = None
T_label = None
t_Label = None

# delay print function
def delay_print(s):
    """Prints a string character by character with a small delay."""
    for c in s:
        sys.stdout.write(c)
        sys.stdout.flush()
        time.sleep(0.010)
    sys.stdout.write('\n')

# checks value is a compatible float
def check_value_float(value):
    try:
        value = float(value)
        return True 
    except ValueError:
        return False 
    
# calculating new value for atmospheric temperature 
def atmosphere_temp(E, y_displacement):
    factor = 9.8/1000
    E = E - (y_displacement * factor)
    return E 

# value finding for the value of Force drag
def calculations_Fdrag(density, velocity, area, Cd):
    p = density  # FIXED: was 'desnity' (typo)
    v = velocity
    A = area
    # Prevent overflow by capping velocity
    if abs(v) > 10000:
        v = 10000 if v > 0 else -10000
    Fdrag = (p * (v**2) * A * Cd) / 2
    # Cap Fdrag to prevent overflow
    if Fdrag > 1e10:
        Fdrag = 1e10
    return Fdrag

# value finding for the value of calculations of deceleration and new value of v in both x and y axis 
def calculaion_rate_deceleration(m, x_velocity, y_velocity, velocity, g, time_period, Fdrag):
    # value for deceleration in the x axis 
    vx = x_velocity
    vy = y_velocity
    v = velocity 
    
    # Prevent division by zero
    if abs(v) < 0.001:
        v = 0.001
    
    # Prevent division by zero mass
    if m < 0.001:
        m = 0.001
    
    # FIXED: Correct physics formulas
    # Drag deceleration components
    drag_decel_x = (-Fdrag * (vx / v)) / m
    drag_decel_y = (-Fdrag * (vy / v)) / m
    
    # Total deceleration (drag + gravity in y direction)
    deceleration_x = drag_decel_x
    deceleration_y = drag_decel_y - g  # FIXED: gravity is acceleration, not force
    
    # Cap accelerations to prevent overflow
    if abs(deceleration_x) > 1000:
        deceleration_x = 1000 if deceleration_x > 0 else -1000
    if abs(deceleration_y) > 1000:
        deceleration_y = 1000 if deceleration_y > 0 else -1000
    
    # Update velocities
    new_vx = vx + (deceleration_x * time_period)
    new_vy = vy + (deceleration_y * time_period)
    
    # Cap velocities to prevent overflow
    if abs(new_vx) > 10000:
        new_vx = 10000 if new_vx > 0 else -10000
    if abs(new_vy) > 10000:
        new_vy = 10000 if new_vy > 0 else -10000
    
    new_v_value = sqrt((new_vx**2) + (new_vy**2))
    
    # Cap total velocity
    if new_v_value > 10000:
        new_v_value = 10000
    
    return new_v_value, new_vx, new_vy  # FIXED: return all values

# calculations for finding the value for the change in temperature of the nose due to drag force
def temperature_intial(Fdrag, v, time_period, m, specific_heat_capacity):
    t = time_period 
    c = specific_heat_capacity
    
    # Prevent division by zero
    if m < 0.001:
        m = 0.001
    if c < 0.001:
        c = 0.001
    
    T = (Fdrag * v * t) / (m * c)
    
    # Cap temperature to prevent overflow
    if T > 1e6:
        T = 1e6
    
    return T

# calculations for temperature transfer to external environment and final temperature of object 
def temperture_actual(T, e, time_period, E, sigma, Fdrag, v, m, specific_heat_capacity):
    t = time_period
    Te = E
    s = sigma 
    c = specific_heat_capacity
    
    # Prevent division by zero
    if m < 0.001:
        m = 0.001
    if c < 0.001:
        c = 0.001
    
    # Prevent overflow in temperature difference calculation
    temp_diff = T - Te
    if abs(temp_diff) > 1000:
        temp_diff = 1000 if temp_diff > 0 else -1000
    
    # Cap the fourth power calculation
    temp_diff_4 = temp_diff**4
    if abs(temp_diff_4) > 1e12:
        temp_diff_4 = 1e12 if temp_diff_4 > 0 else -1e12
    
    Ta = ((Fdrag * v * t) - (e * s * temp_diff_4 * t)) / (m * c)
    
    # Cap final temperature
    if abs(Ta) > 1e6:
        Ta = 1e6 if Ta > 0 else -1e6
    
    return Ta

# subroutine to calculate velocity in x 
def calculations_velocity_x(velocity, a):
    v = velocity 
    r_a = np.radians(a)   
    v_x = velocity * math.cos(r_a)
    return v_x

# subroutine to calculate velocity in y  
def calculations_velocity_y(velocity, a):  # FIXED: was using 'angle' instead of 'a'
    v = velocity 
    r_a = np.radians(a) 
    v_y = velocity * math.sin(r_a)
    return v_y

# checks when the value for Y_max is achieved 
def check_Y_max(velocity_y, time_period, current_displacement_y):
    V_y = velocity_y
    t = time_period
    D_y = current_displacement_y
    if V_y <= 0:
        Y_max = D_y
        return True
    else:
        return False

# checks an iteration loop to find the range 
def find_range(displacement_y, displacement_x):
    d_y = displacement_y
    d_x = displacement_x
    if d_y <= 0:
        range_object = d_x
        return range_object
    else:
        return d_x

# error display panel
def error_numeric():
    messagebox.showerror("Incorrect Error", "Error: please enter valid number")

# collecting variables from gui 
# collecting velocity 
def my_function_v():
    global entry_velocity
    value = entry_velocity.get()
    correct = check_value_float(value)
    if correct == True:
        v = value
        return float(v)  # FIXED: convert to float
    else:
        return False 

# collecting mass value 
def my_function_m():
    global entry_mass
    value = entry_mass.get()
    correct = check_value_float(value)
    if correct == True:
        m = value
        return float(m)  # FIXED: convert to float
    else:
        return False 

# collecting angle value 
def my_function_a():
    global entry_angle
    value = entry_angle.get()
    correct = check_value_float(value)
    if correct == True:
        a = value
        return float(a)  # FIXED: convert to float
    else:
        return False
    
# collecting emissivity value 
def my_function_e():
    global entry_emistivity
    value = entry_emistivity.get()
    correct = check_value_float(value)
    if correct == True:
        e = value
        return float(e)  # FIXED: convert to float
    else:
        return False 
    
# collecting specific heat value 
def my_function_c():
    global entry_specific_heat_capacity
    value = entry_specific_heat_capacity.get()
    correct = check_value_float(value)
    if correct == True:
        c = value
        return float(c)  # FIXED: convert to float
    else:
        return False 
    
# collecting area value 
def my_function_A():
    global entry_area
    value = entry_area.get()
    correct = check_value_float(value)
    if correct == True:
        A = value
        return float(A)  # FIXED: convert to float
    else:
        return False 

# collecting density value 
def my_function_p():
    global entry_density
    value = entry_density.get()
    correct = check_value_float(value)
    if correct == True:
        p = value
        return float(p)  # FIXED: convert to float
    else:
        return False 
    
# collecting drag coefficient 
def my_function_Cd():
    global entry_drag_coeffecient
    value = entry_drag_coeffecient.get()
    correct = check_value_float(value)
    if correct == True:
        Cd = value
        return float(Cd)  # FIXED: convert to float
    else:
        return False   

def my_function_E():
    global entry_enviroment_temperature
    value = entry_enviroment_temperature.get()
    correct = check_value_float(value)
    if correct == True:
        E = value
        return float(E)  # FIXED: convert to float
    else:
        return False   
    
# update subroutine for gui  
# change subroutine for the velocity value
def change_v(v, velocity_Label):
    velocity_Label.config(text="velocity " + str(v) + " ms^-1")  # FIXED: added str()

# change subroutine for the mass value
def change_m(m, mass_Label):
    mass_Label.config(text="mass " + str(m) + " kg")  # FIXED: added str()

# change subroutine for the angle value
def change_a(a, angle_label):
    angle_label.config(text="angle " + str(a) + " degrees")  # FIXED: added str()

# change subroutine for the emissivity 
def change_e(e, emistivity_Label):
    emistivity_Label.config(text="emissivity " + str(e))  # FIXED: added str()

# change specific heat capacity 
def change_c(c, specific_heat_capacity_Label):
    specific_heat_capacity_Label.config(text="specific heat capacity " + str(c) + "")  # FIXED: added str()

# change area
def change_A(A, aera_Label):
    aera_Label.config(text="area " + str(A) + " m^2")  # FIXED: added str(), changed m^3 to m^2

# change density 
def change_p(p, density_Label):
    density_Label.config(text="density " + str(p) + " kg/m^3")  # FIXED: added str(), fixed unit

# change Fdrag value
def change_Fdrag(Fdrag, Fdrag_Label):
    Fdrag_Label.config(text="force in drags is " + str(Fdrag) + " N")  # FIXED: removed Label., added str()

# change Temperature value 
def change_T(T, T_label):
    T_label.config(text="Temperature of the object is " + str(T) + " degrees C")  # FIXED: removed label., added str()

# change environment value 
def change_E(E, enviroment_temperature_label):
    enviroment_temperature_label.config(text="temperature of environment is " + str(E) + " K")  # FIXED: removed E_Label =, added str()

# change time value 
def change_t(t, t_Label):
    t_Label.config(text="time period of the object is " + str(t) + " seconds")  # FIXED: added str()

def value_set_up(root):  # FIXED: removed parameters, pass root instead
    global entry_velocity, entry_mass, entry_angle, entry_emistivity
    global entry_specific_heat_capacity, entry_area, entry_density
    global entry_drag_coeffecient, entry_enviroment_temperature
    global velocity_Label, mass_Label, angle_Label, emistivity_Label
    global specific_heat_capacity_Label, aera_Label, density_Label
    global drag_coeffecient, enviroment_temperature
    
    # Setting up labels of inputted variables
    title_label = Label(root, text="Variables")  
    velocity_Label = Label(root, text="velocity: not set")
    mass_Label = Label(root, text="mass: not set")  
    angle_Label = Label(root, text="angle: not set")
    emistivity_Label = Label(root, text="emissivity: not set")
    specific_heat_capacity_Label = Label(root, text="specific heat capacity: not set")
    aera_Label = Label(root, text="area: not set")
    density_Label = Label(root, text="density: not set")
    drag_coeffecient = Label(root, text="drag coefficient: not set")
    enviroment_temperature = Label(root, text="temperature of your environment: not set")  # FIXED: added closing parenthesis
    
    # Displaying labels
    title_label.grid(row=1, column=0)
    velocity_Label.grid(row=2, column=0)
    mass_Label.grid(row=3, column=0)
    angle_Label.grid(row=4, column=0)
    emistivity_Label.grid(row=5, column=0)
    specific_heat_capacity_Label.grid(row=6, column=0)
    aera_Label.grid(row=7, column=0)
    density_Label.grid(row=8, column=0)
    drag_coeffecient.grid(row=9, column=0)
    enviroment_temperature.grid(row=10, column=0)
    
    # Setting up response zones 
    entry_velocity = tk.Entry(root)
    entry_velocity.grid(row=2, column=1)
    entry_mass = tk.Entry(root)
    entry_mass.grid(row=3, column=1)
    entry_angle = tk.Entry(root)
    entry_angle.grid(row=4, column=1)
    entry_emistivity = tk.Entry(root)
    entry_emistivity.grid(row=5, column=1)
    entry_specific_heat_capacity = tk.Entry(root)
    entry_specific_heat_capacity.grid(row=6, column=1)
    entry_area = tk.Entry(root)
    entry_area.grid(row=7, column=1)
    entry_density = tk.Entry(root)
    entry_density.grid(row=8, column=1)
    entry_drag_coeffecient = tk.Entry(root)
    entry_drag_coeffecient.grid(row=9, column=1)
    entry_enviroment_temperature = tk.Entry(root)
    entry_enviroment_temperature.grid(row=10, column=1)
    
    # Setting up buttons for entering values 
    
    # Entry for velocity 
    def velocity_entered():
        result = my_function_v()
        if result == False:
            error_numeric()
        else:
            change_v(result, velocity_Label)
    
    velocity_button = tk.Button(root, text="Enter", command=velocity_entered)
    velocity_button.grid(row=2, column=2)
    
    # Entry for mass
    def mass_entered():
        result = my_function_m()
        if result == False:
            error_numeric()
        else:
            change_m(result, mass_Label)
    
    mass_button = tk.Button(root, text="Enter", command=mass_entered)
    mass_button.grid(row=3, column=2)

    # Entry for angle 
    def angle_entered():
        result = my_function_a()
        if result == False:
            error_numeric()
        else:
            change_a(result, angle_Label)
    
    angle_button = tk.Button(root, text="Enter", command=angle_entered)
    angle_button.grid(row=4, column=2)
    
    # Entry for emissivity    
    def emissivity_entered():
        result = my_function_e()
        if result == False:
            error_numeric()
        else:
            change_e(result, emistivity_Label)
    
    emistivity_button = tk.Button(root, text="Enter", command=emissivity_entered)
    emistivity_button.grid(row=5, column=2)
    
    # Entry for specific heat capacity 
    def specific_heat_entered():
        result = my_function_c()
        if result == False:
            error_numeric()
        else:
            change_c(result, specific_heat_capacity_Label)
    
    specific_heat_button = tk.Button(root, text="Enter", command=specific_heat_entered)
    specific_heat_button.grid(row=6, column=2)
        
    # Entry for area 
    def area_entered():
        result = my_function_A()
        if result == False:
            error_numeric()
        else:
            change_A(result, aera_Label)
    
    area_button = tk.Button(root, text="Enter", command=area_entered)
    area_button.grid(row=7, column=2)
        
    # Entry for density
    def density_entered():
        result = my_function_p()
        if result == False:
            error_numeric()
        else:
            change_p(result, density_Label)
    
    density_button = tk.Button(root, text="Enter", command=density_entered)
    density_button.grid(row=8, column=2)
        
    # Entry for drag coefficient
    def drag_coeff_entered():
        result = my_function_Cd()
        if result == False:
            error_numeric()
        else:
            drag_coeffecient.config(text="drag coefficient " + str(result))
    
    drag_coeffecient_button = tk.Button(root, text="Enter", command=drag_coeff_entered)
    drag_coeffecient_button.grid(row=9, column=2)
    
    # Entry for environmental temperature 
    def env_temp_entered():
        result = my_function_E()
        if result == False:
            error_numeric()
        else:
            change_E(result, enviroment_temperature)
    
    enviroment_temperature_button = tk.Button(root, text="Enter", command=env_temp_entered)
    enviroment_temperature_button.grid(row=10, column=2)

def calculated_value_set_up(root):  # FIXED: removed Fdrag, T, t parameters
    global Fdrag_Label, T_label, t_Label
    
    # Setting up labels 
    Fdrag_Label = Label(root, text="force in drags is: not calculated")
    T_label = Label(root, text="Temperature of the object is: not calculated")
    t_Label = Label(root, text="time period of the object is: not calculated")
    
    # Displaying labels 
    Fdrag_Label.grid(row=11, column=0)
    T_label.grid(row=12, column=0)
    t_Label.grid(row=13, column=0)

# Iteration subroutine
def calculation_iteration_subroutine(a, m, e, c, p, v, A, Cd, E, 
                                    ax1, ax2, ax3, ax4, ax5, ax6, 
                                    canvas1, canvas2, canvas3, canvas4, canvas5, canvas6, root):
    
    # Validate input parameters to prevent overflow
    if v <= 0 or v > 10000:
        messagebox.showerror("Invalid Input", "Velocity must be between 0 and 10000 m/s")
        return
    if m <= 0 or m > 1e6:
        messagebox.showerror("Invalid Input", "Mass must be between 0 and 1,000,000 kg")
        return
    if a < 0 or a > 90:
        messagebox.showerror("Invalid Input", "Angle must be between 0 and 90 degrees")
        return
    if p <= 0 or p > 10:
        messagebox.showerror("Invalid Input", "Density must be between 0 and 10 kg/m³")
        return
    
    # Lists to store data for plotting
    x_positions = [0]
    y_positions = [0]
    time_values = [0]
    drag_values = [0]
    thermal_values = [0]
    velocity_values = [v]
    x_velocity_values = []
    y_velocity_values = []
    
    # FIXED: Calculate initial velocity components ONCE at the start
    velocity = v
    x_velocity = calculations_velocity_x(v, a)
    y_velocity = calculations_velocity_y(v, a)
    
    x_velocity_values.append(x_velocity)
    y_velocity_values.append(y_velocity)
    
    x_displacement = 0
    y_displacement = 0
    
    t = 0
    time_period = 0.01
    g = 9.81
    sigma = 5.67e-8
    
    max_iterations = 100000  # Prevent infinite loops
    iteration_count = 0
    
    # Iteration loop with error handling
    try:
        # FIXED: Check y_displacement > -1 to allow landing, not >= 0
        while y_displacement > -1 and iteration_count < max_iterations:
            iteration_count += 1
            
            # FIXED: Don't recalculate velocities from angle - use the updated values from deceleration
            # The velocities are updated by the deceleration function below
            
            # Calculate displacement using current velocities
            x_displacement_new = x_velocity * time_period
            y_displacement_new = y_velocity * time_period
            
            # Check for overflow in displacement
            if abs(x_displacement_new) > 1e6 or abs(y_displacement_new) > 1e6:
                messagebox.showwarning("Simulation Stopped", "Displacement values too large - simulation stopped")
                break
            
            x_displacement = x_displacement + x_displacement_new
            y_displacement = y_displacement + y_displacement_new
            
            # FIXED: Break if we've landed (crossed y=0 going down)
            if y_displacement < 0:
                break
            
            # Calculating Fdrag using current total velocity
            Fdrag = calculations_Fdrag(p, velocity, A, Cd) 
            
            # Check for NaN or Inf
            if not np.isfinite(Fdrag):
                messagebox.showerror("Calculation Error", "Invalid drag force calculated")
                break
            
            # Calculating the initial temperature increase of the material
            T = temperature_intial(Fdrag, velocity, time_period, m, c)
            
            # Calculating actual temperature
            Ta = temperture_actual(T, e, time_period, E, sigma, Fdrag, velocity, m, c)
            
            # Check for NaN or Inf in temperature
            if not np.isfinite(Ta):
                Ta = 0
            
            # FIXED: Calculate new velocities based on drag and gravity
            new_velocity, new_x_velocity, new_y_velocity = calculaion_rate_deceleration(
                m, x_velocity, y_velocity, velocity, g, time_period, Fdrag)
            
            # Check for NaN or Inf in velocities
            if not np.isfinite(new_velocity):
                messagebox.showerror("Calculation Error", "Invalid velocity calculated")
                break
            
            # Update velocities for next iteration
            velocity = new_velocity
            x_velocity = new_x_velocity
            y_velocity = new_y_velocity
            
            # Append current values
            x_positions.append(x_displacement)
            y_positions.append(y_displacement)
            time_values.append(t)
            drag_values.append(Fdrag)
            thermal_values.append(Ta)
            velocity_values.append(velocity)
            x_velocity_values.append(x_velocity)
            y_velocity_values.append(y_velocity)
            
            # Update plots every 10 iterations for performance
            if len(time_values) % 10 == 0:
                # Plot 1: Trajectory
                ax1.clear()
                ax1.plot(x_positions, y_positions, 'b-')
                ax1.set_title("Displacement over time")
                ax1.set_xlabel("Range (m)")
                ax1.set_ylabel("Height (m)")
                ax1.grid(True)
                
                # Plot 2: Thermal energy
                ax2.clear()
                ax2.plot(time_values, thermal_values, 'r-')
                ax2.set_title("Thermal energy transferred")
                ax2.set_xlabel("Time (s)")
                ax2.set_ylabel("Temperature change (K)")
                ax2.grid(True)
                
                # Plot 3: Drag force
                ax3.clear()
                ax3.plot(time_values, drag_values, 'g-')
                ax3.set_title("Air resistance/drag over time")
                ax3.set_xlabel("Time (s)")
                ax3.set_ylabel("Drag force (N)")
                ax3.grid(True)
                
                # Plot 4: Y velocity
                ax4.clear()
                ax4.plot(time_values, y_velocity_values, 'purple')
                ax4.set_title("Velocity in y against time")
                ax4.set_xlabel("Time (s)")
                ax4.set_ylabel("Velocity in Y (m/s)")
                ax4.grid(True)
                
                # Plot 5: X velocity
                ax5.clear()
                ax5.plot(time_values, x_velocity_values, 'orange')
                ax5.set_title("Velocity in x against time")
                ax5.set_xlabel("Time (s)")
                ax5.set_ylabel("Velocity in X (m/s)")
                ax5.grid(True)
                
                # Plot 6: Total velocity
                ax6.clear()
                ax6.plot(time_values, velocity_values, 'cyan')
                ax6.set_title("Velocity against time")
                ax6.set_xlabel("Time (s)")
                ax6.set_ylabel("Velocity (m/s)")
                ax6.grid(True)
                
                # Refresh canvases
                canvas1.draw()
                canvas2.draw()
                canvas3.draw()
                canvas4.draw()
                canvas5.draw()
                canvas6.draw()
                root.update()
                
                # Update labels
                change_Fdrag(Fdrag, Fdrag_Label)
                change_t(t, t_Label)
                change_T(Ta, T_label)
            
            t += time_period
            
            # Safety check to prevent excessively long simulations
            if t > 1000:
                messagebox.showinfo("Simulation Complete", "Maximum simulation time (1000s) reached")
                break
        
        if iteration_count >= max_iterations:
            messagebox.showwarning("Simulation Stopped", "Maximum iterations reached")
        
        # Final plot update
        ax1.clear()
        ax1.plot(x_positions, y_positions, 'b-')
        ax1.set_title("Displacement over time")
        ax1.set_xlabel("Range (m)")
        ax1.set_ylabel("Height (m)")
        ax1.grid(True)
        canvas1.draw()
        
        # FIXED: Better completion message
        max_height = max(y_positions) if y_positions else 0
        messagebox.showinfo("Simulation Complete", 
                          f"Simulation finished!\nRange: {x_displacement:.2f} m\nMax Height: {max_height:.2f} m\nTime: {t:.2f} s")
    
    except Exception as ex:
        messagebox.showerror("Calculation Error", f"An error occurred during calculation:\n{str(ex)}")
        return

# Check values of entry before calculation 
def Check_Values(a, m, e, c, p, v, A, Cd, E):
    start = False
    total_correct = 0
    if a != 0:
        total_correct = total_correct + 1
    if m != 0:
        total_correct = total_correct + 1
    if e != 0:
        total_correct = total_correct + 1
    if c != 0:
        total_correct = total_correct + 1
    if p != 0:
        total_correct = total_correct + 1
    if v != 0:
        total_correct = total_correct + 1
    if A != 0:
        total_correct = total_correct + 1
    if Cd != 0:
        total_correct = total_correct + 1
    if E != 0:
        total_correct = total_correct + 1
    if total_correct == 9:
        start = True
    else:
        start = False 
    return start 

# Main execution 
def main():
    # Setting up backdrop 
    root = tk.Tk()
    root.geometry("1600x900")
    root.title("Projectile Motion Simulator")
    root.configure(bg='black')
    
    # Font library 
    font1 = {'family': 'sans-serif', 'color': 'black', 'size': 15}
    font2 = {'family': 'sans-serif', 'color': 'black', 'size': 10}  

    # Create first figure Displacement (top-left)
    fig1 = Figure(figsize=(4, 2.5))
    fig1.patch.set_facecolor('grey')
    ax1 = fig1.add_subplot(111)
    ax1.set_title("Displacement over time")
    ax1.set_xlabel("Range (m)")
    ax1.set_ylabel("Height (m)")
    ax1.plot(0, 0)
    ax1.grid(True)
    
    canvas1 = FigureCanvasTkAgg(fig1, master=root)
    canvas1.draw()
    canvas1.get_tk_widget().grid(row=0, column=3, rowspan=5, columnspan=2, sticky='nsew', padx=5, pady=5)
    
    # Create second figure Thermal energy (top-middle)
    fig2 = Figure(figsize=(4, 2.5))
    fig2.patch.set_facecolor('grey')
    ax2 = fig2.add_subplot(111)
    ax2.set_title("Thermal energy transferred")
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Temperature change (K)")
    ax2.plot(0, 0)
    ax2.grid(True)
    
    canvas2 = FigureCanvasTkAgg(fig2, master=root)
    canvas2.draw()
    canvas2.get_tk_widget().grid(row=0, column=5, rowspan=5, columnspan=2, sticky='nsew', padx=5, pady=5)
    
    # Create third figure air resistance (top-right)
    fig3 = Figure(figsize=(4, 2.5))
    fig3.patch.set_facecolor('grey')
    ax3 = fig3.add_subplot(111)
    ax3.set_title("Air resistance/drag over time")
    ax3.set_xlabel("Time (s)")
    ax3.set_ylabel("Drag force (N)")
    ax3.plot(0, 0)
    ax3.grid(True)
    
    canvas3 = FigureCanvasTkAgg(fig3, master=root)
    canvas3.draw()
    canvas3.get_tk_widget().grid(row=0, column=7, rowspan=5, columnspan=2, sticky='nsew', padx=5, pady=5)
    
    # Create fourth figure (bottom-left)
    fig4 = Figure(figsize=(4, 2.5))
    fig4.patch.set_facecolor('grey')
    ax4 = fig4.add_subplot(111)
    ax4.set_title("Velocity in y against time")
    ax4.set_xlabel("Time (s)")
    ax4.set_ylabel("Velocity (m/s)")
    ax4.plot(0, 0)
    ax4.grid(True)
    
    canvas4 = FigureCanvasTkAgg(fig4, master=root)
    canvas4.draw()
    canvas4.get_tk_widget().grid(row=5, column=3, rowspan=5, columnspan=2, sticky='nsew', padx=5, pady=5)
    
    # Create fifth figure (bottom-middle)
    fig5 = Figure(figsize=(4, 2.5))
    fig5.patch.set_facecolor('grey')
    ax5 = fig5.add_subplot(111)
    ax5.set_title("Velocity in x against time")
    ax5.set_xlabel("Time (s)")
    ax5.set_ylabel("Velocity (m/s)")
    ax5.plot(0, 0)
    ax5.grid(True)
    
    canvas5 = FigureCanvasTkAgg(fig5, master=root)
    canvas5.draw()
    canvas5.get_tk_widget().grid(row=5, column=5, rowspan=5, columnspan=2, sticky='nsew', padx=5, pady=5)
    
    # Create sixth figure (bottom-right)
    fig6 = Figure(figsize=(4, 2.5))
    fig6.patch.set_facecolor('grey')
    ax6 = fig6.add_subplot(111)
    ax6.set_title("Velocity against time")
    ax6.set_xlabel("Time (s)")
    ax6.set_ylabel("Velocity (m/s)")
    ax6.plot(0, 0)
    ax6.grid(True)
    
    canvas6 = FigureCanvasTkAgg(fig6, master=root) 
    canvas6.draw()
    canvas6.get_tk_widget().grid(row=5, column=7, rowspan=5, columnspan=2, sticky='nsew', padx=5, pady=5)
    
    # Configure grid weights to make graphs expand properly
    for i in range(3, 9):
        root.columnconfigure(i, weight=1)
    for i in range(10):
        root.rowconfigure(i, weight=1)
    
    # Setting up the labels for the variables on the screen 
    value_set_up(root)
    calculated_value_set_up(root)
    
    # Calculate button
    def calculate():
        # Get current values from labels
        v = my_function_v() or 0
        m = my_function_m() or 0
        a = my_function_a() or 0
        e = my_function_e() or 0
        c = my_function_c() or 0
        A = my_function_A() or 0
        p = my_function_p() or 0
        Cd = my_function_Cd() or 0
        E = my_function_E() or 0
        
        start_possible = Check_Values(a, m, e, c, p, v, A, Cd, E)
        if start_possible == True:
            calculation_iteration_subroutine(a, m, e, c, p, v, A, Cd, E,
                                            ax1, ax2, ax3, ax4, ax5, ax6, 
                                            canvas1, canvas2, canvas3, canvas4, canvas5, canvas6, root)
        else:
            messagebox.showerror("Missing Values", "Please enter all values before calculating")
    
    calculate_button = tk.Button(root, text="Calculate", command=calculate, bg='green', fg='white', font=('Arial', 14))
    calculate_button.grid(row=14, column=1, pady=10)
    
    root.mainloop()

if __name__ == "__main__":
    main()
