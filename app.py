from flask import Flask, render_template, request, Response
from math import sqrt, pow, radians, degrees
import numpy as np
from sympy import symbols, Eq, solve, sin, cos, tan, atan, nsolve
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import io
import os
from scipy.optimize import fsolve
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"

app = Flask(__name__)

#-------------------------------------------------------------------
def isentropic_calc(var1, value1, gamma, R):
    results = {}
    M = symbols('M')    
    if var1 == "M":
        M = value1
    elif var1 == "Mstar":
        M = float(nsolve(Eq(((((M**2) * (gamma + 1)) / (2 + (gamma - 1) * (M**2)))**0.5), value1), M, 1))
    elif var1 == "P0_P":
        M = float(nsolve(Eq((1 + (gamma - 1) * 0.5 * (M ** 2))**(gamma / (gamma - 1)), value1), M, 1))
    elif var1 == "T0_T":
        M = float(nsolve(Eq((1 + (gamma - 1) * 0.5 * (M ** 2)), value1), M, 1))
    elif var1 == "rho0_rho":
        M = float(nsolve(Eq((1 + (gamma - 1) * 0.5 * (M ** 2))**(1 / (gamma - 1)), value1), M, 1))
    
    P0_P = pow((1 + (gamma - 1) * 0.5 * (M ** 2)), (gamma / (gamma - 1)))
    rho0_rho = pow((1 + (gamma - 1) * 0.5 * (M ** 2)), (1 / (gamma - 1)))
    T0_T = (1 + (gamma - 1) * 0.5 * (M ** 2))
    Mstar = sqrt(((M**2) * (gamma + 1)) / (2 + (gamma - 1) * (M**2)))
    Pstar_P = pow((2 / (gamma + 1)), (gamma / (gamma - 1)))
    rhostar_rho = pow((2 / (gamma + 1)), (1 / (gamma - 1)))
    Tstar_T = (2 / (gamma + 1))
    astar_a0 = sqrt(2 / (gamma + 1))
    A_Astar = sqrt(((M**2) * (gamma + 1)) / (2 + (gamma - 1) * (M**2)))    
    if M is not None:
        results = {
            "M" : f"{M:.4f}",
            "P₀/P" : f"{P0_P:.4f}",
            "ρ₀/ρ" : f"{rho0_rho:.4f}",
            "T₀/T" : f"{T0_T:.4f}",
            "M*" : f"{Mstar:.4f}",
            "P*/P" : f"{Pstar_P:.4f}",
            "ρ*/ρ" : f"{rhostar_rho:.4f}",
            "T*/T" : f"{Tstar_T:.4f}",
            "a*/a₀" : f"{astar_a0:.4f}",
            "A/A*" : f"{A_Astar:.4f}"
        }    
    return results, float(M)

#-------------------------------------------------------------------
def normal_calc(var1, value1, gamma, R):
    M1 = symbols('M1')
    equations = {
        "M1": value1,
        "M2": Eq(value1**2, (1 + ((gamma - 1) / 2) * M1**2) / (gamma * M1**2 - (gamma - 1) / 2)),
        "P2_P1": Eq(value1, 1 + (2 * gamma * (M1**2 - 1)) / (gamma + 1)),
        "rho2_rho1": Eq(value1, ((gamma +1) * (M1**2)) / (2 + (gamma - 1) * (M1**2))),
        "T2_T1": Eq(value1, (1 + (2 * gamma * (M1**2 - 1)) / (gamma + 1)) / (( (gamma + 1) * M1**2) / (2 + (gamma - 1) * M1**2)))
    }
    if var1 not in equations:
        return {"error": "Invalid var"}
    if isinstance(equations[var1], Eq):
        try:
            M1 = float(nsolve(equations[var1], M1, 1))
        except:
            return {"error": "Could not solve for M1"}
    else:
        M1 =float(value1)

    M2 = ((1 + ((gamma - 1) / 2) * M1**2) / (gamma * M1**2 - (gamma - 1) / 2))**0.5
    P2_P1 = 1 + (2 * gamma * (M1**2 - 1)) / (gamma + 1)
    rho2_rho1 = ((gamma +1) * (M1**2)) / (2 + (gamma - 1) * (M1**2))
    T2_T1 = P2_P1 / rho2_rho1
    P02_P01 = ((gamma + 1) * M1**2 / ((gamma - 1) * M1**2 + 2))**(gamma / (gamma - 1)) * ((gamma + 1) / (2 * gamma * M1**2 - (gamma - 1)))**(1 / (gamma - 1))    
    results = {
        "M₁": f"{M1:.4f}",
        "M₂": f"{M2:.4f}",
        "P₂/P₁": f"{P2_P1:.4f}",
        "ρ₂/ρ₁": f"{rho2_rho1:.4f}",
        "T₂/T₁": f"{T2_T1:.4f}",
        "P₀₂/P₀₁": f"{P02_P01:.4f}",
        "ρ₀₂/ρ₀₁": f"{P02_P01:.4f}",
        "T₀₂/T₀₁": f"{1:.4f}"
    }    
    return results

#-------------------------------------------------------------------
def oblique_calc(var1, value1, var2, value2, gamma, R):
    results = {}
    M1 = symbols('M1')
    theta = symbols('theta')
    beta = symbols('beta')
    if var1 == "M1":
        M1 = value1
        if var2 == "theta":
            theta = radians(value2) 
            beta_guess = radians(20) 
            beta_sol = nsolve(Eq(theta, atan((2 / tan(beta)) * (((M1**2 * sin(beta)**2 - 1) / (M1**2 * (gamma + cos(2 * beta)) + 2))))), beta, beta_guess)
            if beta_sol is None:
                return {"error": "No valid oblique shock solution found"}
            beta = float(beta_sol)
        elif var2 == "beta":
            beta = radians(value2)  
            theta_sol = solve(Eq(tan(theta), (2 / tan(beta)) * 
                        (((M1**2 * sin(beta)**2 - 1) / 
                        (M1**2 * (gamma + cos(2 * beta)) + 2)))), theta)
            theta = float(theta_sol[0].evalf())
        Mn1 = M1 * sin(beta)
        Mn2 = sqrt((1 + ((gamma - 1) / 2) * Mn1**2) / (gamma * Mn1**2 - (gamma - 1) / 2))
        M2 = Mn2 / sin(beta - theta) 
        
    if M1 is not None and theta is not None and beta is not None: 
        results["M₁"] = f"{M1:.4f}"
        results["M₂"] = f"{M2:.4f}"
        results["Mₙ₁"] = f"{Mn1:.4f}"
        results["Mₙ₂"] = f"{Mn2:.4f}"
        results["θ"] = f"{degrees(theta):.4f}"
        results["β"] = f"{degrees(beta):.4f}"
        results["P₂/P₁"] = f"{1 + (2 * gamma * ((Mn1**2) - 1) / (gamma + 1)):.4f}"
        results["ρ₂/ρ₁"] = f"{(gamma+1) * Mn1**2 / (2 + (gamma-1) * Mn1**2):.4f}"
        results["T₂/T₁"] = f"{((1 + (2 * gamma * ((Mn1**2) - 1) / (gamma + 1))) / ((gamma+1) * Mn1**2 / (2 + (gamma-1) * Mn1**2))):.4f}"
    return results

#-------------------------------------------------------------------
@app.route("/plot/isentropic/<float:M>")
def plot_isentropic(M):
    gamma = 1.4
    M_min = max(0.1, M - 3)
    M_max = M + 3
    M_vals = np.linspace(M_min, M_max, 100)
    P0_P = (1 + (gamma - 1) * 0.5 * (M_vals ** 2)) ** (gamma / (gamma - 1))
    T0_T = (1 + (gamma - 1) * 0.5 * (M_vals ** 2))
    rho0_rho = (1 + (gamma - 1) * 0.5 * (M_vals ** 2)) ** (1 / (gamma - 1))

    plt.figure(figsize=(5, 4))
    plt.plot(M_vals, P0_P, label="P/P₀")
    plt.plot(M_vals, rho0_rho, label="ρ/ρ₀")
    plt.plot(M_vals, T0_T, label="T/T₀")
    plt.xlabel("Mach Number (M)")
    plt.ylabel("Ratios")
    plt.legend()
    plt.grid()
    buf = io.BytesIO()
    plt.savefig(buf, format="png")
    buf.seek(0)
    plt.close()
    return Response(buf.getvalue(), mimetype="image/png")

#-------------------------------------------------------------------
@app.route("/plot/normal")
def plot_normal():
    gamma = 1.4
    M1 = np.linspace(0.5, 10, 100) 
    M2 = np.sqrt(np.divide(1 + ((gamma - 1) / 2) * M1**2, gamma * M1**2 - (gamma - 1) / 2))
    P2_P1 = 1 + (2 * gamma * (M1**2 - 1)) / (gamma + 1)
    rho2_rho1 = ((gamma + 1) * M1**2) / (2 + (gamma - 1) * M1**2)
    T2_T1 = np.divide(P2_P1, rho2_rho1) 
    P02_P01 = np.power(((gamma + 1) * M1**2) / ((gamma - 1) * M1**2 + 2), gamma / (gamma - 1)) * np.power((gamma + 1) / (2 * gamma * M1**2 - (gamma - 1)), 1 / (gamma - 1))

    plt.figure(figsize=(5, 4))
    plt.plot(M1, M2, label="M₂ (Downstream Mach)", linestyle='--')
    plt.plot(M1, P2_P1/100, label="P₂/P₁ (Pressure Ratio)", linestyle='-')
    plt.plot(M1, rho2_rho1/10, label="ρ₂/ρ₁ (Density Ratio)", linestyle='-')
    plt.plot(M1, T2_T1/10, label="T₂/T₁ (Temperature Ratio)", linestyle='-')
    plt.plot(M1, P02_P01, label="P₀₂/P₀₁ (Total Pressure Ratio)", linestyle='-')
    plt.xlabel("Upstream Mach Number (M₁)")
    plt.ylabel("Ratios")
    plt.legend()
    plt.grid()
    buf = io.BytesIO()
    plt.savefig(buf, format="png")
    buf.seek(0)
    plt.close()
    return Response(buf.getvalue(), mimetype="image/png")

#-------------------------------------------------------------------
def solve_beta(M, theta, gamma=1.4):
    theta_rad = np.radians(theta)
    def beta_eq(beta):
        return np.tan(theta_rad) - (2 / np.tan(beta)) * (
            ((M ** 2 * np.sin(beta) ** 2 - 1) / (M ** 2 * (gamma + np.cos(2 * beta)) + 2))
        )
    weak_beta_guess = np.radians(20)
    strong_beta_guess = np.radians(80)
    beta_solutions = []

    for beta_guess in [weak_beta_guess, strong_beta_guess]:
        try:
            beta_sol = fsolve(beta_eq, beta_guess)[0] 
            if beta_sol > 0:
                beta_solutions.append(np.degrees(beta_sol)) 
        except:
            continue 
    return beta_solutions  

#-------------------------------------------------------------------
@app.route("/plot/oblique")
def plot_oblique():
    gamma = 1.4 
    M1, M2 = 2, 3
    theta_deg = np.linspace(0, 40, 100)

    plt.figure(figsize=(5, 4))
    all_data_M1 = [] 
    all_data_M2 = [] 
    for M, data_store in [(M1, all_data_M1), (M2, all_data_M2)]:
        for th in theta_deg:
            beta_values = solve_beta(M, th, gamma)
            for beta in beta_values:
                data_store.append((th, beta)) 
    all_data_M1 = np.array(all_data_M1)
    all_data_M2 = np.array(all_data_M2)
    if all_data_M1.size > 0:
        plt.scatter(all_data_M1[:, 0], all_data_M1[:, 1], color="lightblue", label=f"M = {M1}", s=10, alpha=0.7)
    if all_data_M2.size > 0:
        plt.scatter(all_data_M2[:, 0], all_data_M2[:, 1], color="orange", label=f"M = {M2}", s=10, alpha=0.7)

    plt.xlabel("Deflection Angle (θ)°")
    plt.ylabel("Shock Angle (β)°")
    plt.legend()
    plt.grid()
    buf = io.BytesIO()
    plt.savefig(buf, format="png")
    buf.seek(0)
    plt.close()
    return Response(buf.getvalue(), mimetype="image/png")

#-------------------------------------------------------------------
@app.route("/", methods=["GET", "POST"])
def home():

    results = None
    plot_url = None
    if request.method == "POST":
        flow_type = request.form["flow_type"]
        gamma = float(request.form.get("gamma_val"))
        R = float(request.form.get("R_val"))
        
        if flow_type == "isentropic":
            var1 = request.form["isen_var1"]
            value1 = float(request.form["isen_val1"])            
            results, computed_M = isentropic_calc(var1, value1, gamma, R)
            plot_url = f"/plot/isentropic/{computed_M}"            
        elif flow_type == "normal":
            var1 = request.form["norm_var1"]
            value1 = float(request.form["norm_val1"])
            results = normal_calc(var1, value1, gamma, R)
            plot_url = "/plot/normal"           
        elif flow_type == "oblique":
            var1 = request.form["obl_var1"]
            value1 = float(request.form["obl_val1"])
            var2 = request.form["obl_var2"]
            value2 = float(request.form["obl_val2"])
            results = oblique_calc(var1, value1, var2, value2, gamma, R)
            plot_url = "/plot/oblique"
        else:
            results = {"error": "Invalid flow type"}
            
    return render_template("index.html", results=results, plot_url=plot_url)

if __name__ == "__main__":
    app.run(debug=True)