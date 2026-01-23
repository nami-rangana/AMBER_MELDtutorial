#!/usr/bin/env python3
import sys
import math

def generate_constant(num_replicas, temperature):
    temps = []
    for i in range(num_replicas):
        alpha = i / (num_replicas - 1) if num_replicas > 1 else 0.0
        temps.append(temperature)
    return temps

def generate_linear(num_replicas, alpha_min, alpha_max, temp_min, temp_max):
    temps = []
    delta_alpha = alpha_max - alpha_min
    delta_temp = temp_max - temp_min
    
    for i in range(num_replicas):
        alpha = i / (num_replicas - 1) if num_replicas > 1 else 0.0
        
        if alpha <= alpha_min:
            temps.append(temp_min)
        elif alpha <= alpha_max:
            frac = (alpha - alpha_min) / delta_alpha
            temp = temp_min + frac * delta_temp
            temps.append(temp)
        else:
            temps.append(temp_max)
    
    return temps

def generate_geometric(num_replicas, alpha_min, alpha_max, temp_min, temp_max):
    temps = []
    delta_alpha = alpha_max - alpha_min
    
    for i in range(num_replicas):
        alpha = i / (num_replicas - 1) if num_replicas > 1 else 0.0
        
        if alpha <= alpha_min:
            temps.append(temp_min)
        elif alpha <= alpha_max:
            frac = (alpha - alpha_min) / delta_alpha
            delta = math.log(temp_max) - math.log(temp_min)
            temp = math.exp(delta * frac + math.log(temp_min))
            temps.append(temp)
        else:
            temps.append(temp_max)
    
    return temps

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: ./generate_temperatures.py <method> <num_replicas> <args...>")
        sys.exit(1)
    
    method = sys.argv[1]
    num_replicas = int(sys.argv[2])
    
    if method == "constant":
        temperature = float(sys.argv[3])
        temps = generate_constant(num_replicas, temperature)
    elif method == "linear":
        alpha_min, alpha_max = float(sys.argv[3]), float(sys.argv[4])
        temp_min, temp_max = float(sys.argv[5]), float(sys.argv[6])
        temps = generate_linear(num_replicas, alpha_min, alpha_max, temp_min, temp_max)
    elif method == "geometric":
        alpha_min, alpha_max = float(sys.argv[3]), float(sys.argv[4])
        temp_min, temp_max = float(sys.argv[5]), float(sys.argv[6])
        temps = generate_geometric(num_replicas, alpha_min, alpha_max, temp_min, temp_max)
    else:
        print(f"Unknown method: {method}")
        sys.exit(1)
    
    with open("temperatures.dat", "w") as f:
        for temp in temps:
            f.write(f"{temp:.1f}\n")
    
    print(f"temperatures.dat generated with {num_replicas} replicas using {method} method")