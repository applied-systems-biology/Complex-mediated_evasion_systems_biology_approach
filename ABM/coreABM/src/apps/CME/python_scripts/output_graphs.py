#  Copyright by Yann Bachelot
#
#  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
#  https://www.leibniz-hki.de/en/applied-systems-biology.html
#  HKI-Center for Systems Biology of Infection
#  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Institute (HKI)
#  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
#  The project code is licensed under BSD2-Clause.
#  See the LICENSE file provided with the code for the full license.

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import os
import sys

path = sys.argv[1] # path to output directory

colors=['rgb(31, 119, 180)', 'rgb(255, 127, 14)',
        'rgb(44, 160, 44)', 'rgb(214, 39, 40)',
        'rgb(148, 103, 189)', 'rgb(140, 86, 75)',
        'rgb(227, 119, 194)', 'rgb(127, 127, 127)',
        'rgb(188, 189, 34)']

colors_name=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'grey', 'yellow']

def rgb_to_rgba(rgb_value, alpha):
    """
    Adds the alpha channel to an RGB Value and returns it as an RGBA Value
    :param rgb_value: Input RGB Value
    :param alpha: Alpha Value to add  in range [0,1]
    :return: RGBA Value
    """
    return f"rgba{rgb_value[3:-1]}, {alpha})"


if "uptake.csv" in os.listdir(path):
    df = pd.read_csv(path+"/uptake.csv", sep=";")
    if not(df.empty):
        if "pop.csv" in os.listdir(path):
            df2 = pd.read_csv(path+"/pop.csv", sep=";")
            if not(df2.empty):
                df2 = df2.iloc[:,:-1]
                total_amp = df2.where(df2["time"] == 0.0).dropna().where(df2["agent"]== "AMP").dropna()["nb"][0]
        else:
            total_amp = 1
        df_f = df.groupby("time", as_index=False)["uptake"].mean()
        std = df.groupby("time", as_index=False)["uptake"].std()["uptake"].tolist()
        df_f["norm_uptake"] = df_f["uptake"]/total_amp *100
        df_f["std"] = std/total_amp
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=df_f['time'], y=df_f["norm_uptake"],
                                 line=dict(color=colors[0], dash='solid'),
                                 name = "uptake"))
        fig.add_trace(go.Scatter(x=df_f['time'], y=df_f["norm_uptake"] + df_f["std"],
                                 line=dict(width=0,  color = rgb_to_rgba(colors[0], 0.0)),
                                 showlegend=False))
        fig.add_trace(go.Scatter(x=df_f['time'], y=df_f["norm_uptake"] - df_f["std"],
                                 line=dict(width=0,  color = rgb_to_rgba(colors[0], 0.0)),
                                 fill='tonexty',
                                 fillcolor= rgb_to_rgba(colors[0], 0.3),
                                 showlegend=False))
        fig.update_layout(title = "Uptake over time", yaxis_title = "Percentage of AMPs uptaken", xaxis_title = "time", font = dict(size = 15))
        with open(path+"/uptake.html", "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

if "uptake_rate.csv" in os.listdir(path):
    df = pd.read_csv(path + "/uptake_rate.csv", sep = ";")
    if not(df.empty):
        df_f = df.groupby("time", as_index=False)["uptake_rate"].mean()
        std = df.groupby("time", as_index=False)["uptake_rate"].std()["uptake_rate"].tolist()
        df_f["std"] = std
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=df_f['time'], y=df_f["uptake_rate"],
                                 line=dict(color=colors[0], dash='solid'),
                                 name = "uptake_rate"))
        fig.add_trace(go.Scatter(x=df_f['time'], y=df_f["uptake_rate"] + df_f["std"],
                                 line=dict(width=0,  color = rgb_to_rgba(colors[0], 0.0)),
                                 showlegend=False))
        fig.add_trace(go.Scatter(x=df_f['time'], y=df_f["uptake_rate"] - df_f["std"],
                                 line=dict(width=0,  color = rgb_to_rgba(colors[0], 0.0)),
                                 fill='tonexty',
                                 fillcolor= rgb_to_rgba(colors[0], 0.3),
                                 showlegend=False))
        fig.update_layout(yaxis_title = "Uptake rate (s-1)", xaxis_title = "time")
        with open(path+"/uptake_rate.html", "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

if "lifetime.csv" in os.listdir(path):
    df = pd.read_csv(path + "/lifetime.csv", sep=";")
    if not(df.empty):
        timestep = df.at[0, "timestep"]
        fig = px.histogram(df["lifetime_complex"]*timestep)
        fig.update_layout(title = "Lifetime of Complex", yaxis_title = "Count", xaxis_title = "time")
        fig2 = px.box(df["lifetime_complex"]*timestep)
        fig2.update_layout(title = "Distribution", yaxis_title = "time",font = dict(size = 15))
        with open(path+"/lifetime_complex.html", "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))
            f.write(fig2.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

if "pop.csv" in os.listdir(path):
    df = pd.read_csv(path + "/pop.csv", sep=";")
    if not(df.empty):
        df = df.iloc[:,:-1]
        df_res = df.groupby(["agent", "time"], as_index=False)["nb"].mean()
        if (len(np.unique(df["id"])) > 1):
            std = df.groupby(["agent", "time"], as_index=False)["nb"].std()["nb"].tolist()
        else:
            std = np.zeros(df_res.shape[0])
        df_res["std"] = std
        fig = go.Figure()
        count = 0
        for i in np.unique(df_res["agent"]):
            df_tmp = df_res.where(df_res["agent"] == i).dropna()
            fig.add_trace(go.Scatter(x=df_tmp['time'], y=df_tmp["nb"],
                                     line=dict(color=colors[count], dash='solid', shape='spline'),
                                     name = i))
            fig.add_trace(go.Scatter(x=df_tmp['time'], y=df_tmp["nb"] + df_tmp["std"],
                                     line=dict(width=0, color = rgb_to_rgba(colors[count], 0.0)),
                                     showlegend=False))
            fig.add_trace(go.Scatter(x=df_tmp['time'], y=df_tmp["nb"] - df_tmp["std"],
                                     line=dict(width=0, color = rgb_to_rgba(colors[count], 0.0)),
                                     fill='tonexty',
                                     fillcolor= rgb_to_rgba(colors[count], 0.3),
                                     showlegend=False))
            count += 1
        fig.update_layout(title = "Evolution of the different agent populations", yaxis_title = "Agents concentration (um-3)", xaxis_title = "time", font = dict(size = 15))
        with open(path+"/population_statistics.html", "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

    df = pd.read_csv(os.path.join(path, "pop.csv"), sep=";")

    if not df.empty:
        df = df.iloc[:, :-1]

        fig = go.Figure()
        unique_ids = df["id"].unique()

        count = 0
        for sim_id in unique_ids:
            df_sim = df[df["id"] == sim_id]

            for agent in df_sim["agent"].unique():
                df_agent = df_sim[df_sim["agent"] == agent]
                df_agent_sorted = df_agent.sort_values("time")

                fig.add_trace(go.Scatter(
                    x=df_agent_sorted['time'],
                    y=df_agent_sorted['nb'],
                    line=dict(color=colors[count % len(colors)], dash='solid', shape='spline'),
                    name=f"{agent} - Sim {sim_id}",
                ))

                count += 1

        fig.update_layout(
            title="Evolution of Agent Populations for Individual Simulations",
            yaxis_title="Agents concentration (um-3)",
            xaxis_title="Time",
            font=dict(size=15)
        )

        with open(os.path.join(path, "population_individual_simulations.html"), "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))
"""
if "agent-statistics.csv" in os.listdir(path):
    df = pd.read_csv(path+"/agent-statistics.csv", sep=";")
    if not(df.empty):
        df = df.iloc[:,:-1]
        df["count"] = np.ones(df.shape[0])/len(np.unique(df["id"]))
        df = df[df["agent"] != "Pathogenic_cell"]
        min_dist = df["dist_to_cell"].min()
        #print("Min dist: "+ str(min_dist))
        tmp_max = df["dist_to_cell"].max()
        #print("Max dist: "+ str(tmp_max))
        cube_size = np.round(2*(tmp_max)/np.sqrt(3))
        #print("cube size: " + str(cube_size))
        df = df[df["dist_to_cell"] < 0.7*cube_size]
        max_dist = df["dist_to_cell"].max()
        #print("Max dist: "+ str(max_dist))

        num_compartments = 22
        step = (np.sqrt(3)*cube_size/2 -min_dist)/num_compartments

        # create a list of the compartments
        compartment = []
        for distance in df["dist_to_cell"]:
            # determine the index of the compartment that this distance belongs to
            compartment.append(int((distance - min_dist) / step))
        df["compartment"] = compartment
        compartment_volume = {}
        for i in range(num_compartments):
            r1 = i * step + min_dist
            v_sphere1 = (4/3) * np.pi * (r1**3)
            r2 = (i + 1) * step + min_dist
            v_sphere2 = (4/3) * np.pi * (r2**3)
            if r2 <= cube_size/2:  # sphere is entirely inside the cube
                compartment_volume[i]= (v_sphere2 - v_sphere1)
            else:
                h2 = r2 - (cube_size/2)
                v_cap2 = ((np.pi * h2**2)/3)*(3*r2-h2)
                v_intersect2 = v_sphere2 - (6*v_cap2)
                h1 = r1 - (cube_size/2)
                v_cap1 = ((np.pi * h1**2)/3)*(3*r1-h1)
                v_intersect1 = v_sphere1 - (6*v_cap1)
                v_f = v_intersect2 - v_intersect1
                compartment_volume[i]= v_f
        #compartment_volume[0] = compartment_volume[0]
        #compartment_volume[num_compartments-1] = compartment_volume[num_compartments-1]/4
        df_test = df.groupby(["time", "agent", "compartment"], as_index=False)["count"].sum()
        #print(df_test.head())
        df_test["norm_count"] = df_test.apply(lambda row: row["count"] / compartment_volume[row["compartment"]], axis=1)
        #print(df_test.head())
        #df_test["norm_count"] = df_test["norm_count"]/1000
        df_test["distance"] = df_test["compartment"]*step + min_dist
        tmp = pd.DataFrame([[0.0, "Complex", 0, 0, 0, min_dist]], columns= df_test.columns)
        tmp2 = pd.DataFrame([[0.0, "Defensive", 0, 0, 0, min_dist]], columns= df_test.columns)
        df_test = pd.concat([df_test, tmp, tmp2], ignore_index=True)
        df_test = df_test.where(df_test["agent"] != "Cell").dropna()
        #print(df_test[df_test["agent"]=="Complex"])
        fig = px.line(df_test, x = "distance", y = "norm_count", animation_frame="time", color="agent", markers = True, line_shape="spline", color_discrete_map={
            "AMP": "blue",
            "Complex": "green",
            "Defensive": "red"
        })
        fig.add_vline(x=min_dist, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")
        fig.add_vrect(x0=0, x1=min_dist, line_width=0, line_dash = "dash", fillcolor="orange", opacity=0.2)
        fig.add_annotation(x = min_dist/2, y = df_test["norm_count"].max()/4, text = "Pathogenic cell", showarrow=False)
        fig.update_layout(xaxis_title = r"$\text{Distance to cell (}\mu m\text{)}$", yaxis_title = r"$\large{\text{Concentration (agents }10^{3}/\mu m^{3})}$", legend_title = "Agents", font = dict(size = 15))
        with open(path+"/spatial_info.html", "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn',  auto_play=False))
if "agent-statistics-multiple_receptors.csv" in os.listdir(path):
    df = pd.read_csv(path+"/agent-statistics-multiple_receptors.csv", sep=";")
    if not(df.empty):
        df = df.iloc[:,:-1]
        df["count"] = np.ones(df.shape[0])/len(np.unique(df["id"]))
        df = df[df["agent"] != "Pathogenic_cell"]
        min_dist = df["dist_to_cell"].min()
        tmp_max = df["dist_to_cell"].max()
        cube_size = np.round(2*(tmp_max)/np.sqrt(3))
        df = df[df["dist_to_cell"] < 0.7*cube_size]
        max_dist = df["dist_to_cell"].max()

        num_compartments = 25
        step = (np.sqrt(3)*cube_size/2 -min_dist)/num_compartments

        compartment = []
        for distance in df["dist_to_cell"]:
            compartment.append(int((distance - min_dist) / step))
        df["compartment"] = compartment
        compartment_volume = {}
        for i in range(num_compartments):
            r1 = i * step + min_dist
            v_sphere1 = (4/3) * np.pi * (r1**3)
            r2 = (i + 1) * step + min_dist
            v_sphere2 = (4/3) * np.pi * (r2**3)
            if r2 <= cube_size/2:
                compartment_volume[i]= (v_sphere2 - v_sphere1)
            else:
                h2 = r2 - (cube_size/2)
                v_cap2 = ((np.pi * h2**2)/3)*(3*r2-h2)
                v_intersect2 = v_sphere2 - (6*v_cap2)
                h1 = r1 - (cube_size/2)
                v_cap1 = ((np.pi * h1**2)/3)*(3*r1-h1)
                v_intersect1 = v_sphere1 - (6*v_cap1)
                v_f = v_intersect2 - v_intersect1
                compartment_volume[i]= v_f

        # Filter out rows with "Complex" agent
        df_complex = df[df["agent"] == "Complex"]

        # Add a new column combining "agent" and "nb_AMP"
        df_complex["complex_type"] = df_complex["agent"] + "_" + df_complex["nb_AMP"].astype(str)

        # Group by time, agent type, and compartment, and sum up counts
        df_test_complex = df_complex.groupby(["time", "complex_type", "compartment"], as_index=False)["count"].sum()

        # Calculate normalized count for each compartment
        df_test_complex["norm_count"] = df_test_complex.apply(lambda row: row["count"] / compartment_volume[row["compartment"]], axis=1)

        # Convert "norm_count" to thousandths
        df_test_complex["norm_count"] = df_test_complex["norm_count"] / 1000

        # Create a new dataframe to include all complex types
        df_combined_complex = pd.DataFrame(columns=df_test_complex.columns)

        # Iterate over unique complex types and append to the combined dataframe
        for complex_type in df_test_complex["complex_type"].unique():
            df_tmp = df_test_complex[df_test_complex["complex_type"] == complex_type]
            df_tmp["agent"] = complex_type  # Rename "complex_type" to "agent"
            df_combined_complex = pd.concat([df_combined_complex, df_tmp])

        # Recalculate the "distance" column for df_combined_complex
        df_combined_complex["distance"] = df_combined_complex["compartment"] * step + min_dist

        # Add "AMP" and "Defensive" data to the combined dataframe
        df_AMP = df[df["agent"] == "AMP"]
        df_defensive = df[df["agent"] == "Defensive"]

        df_combined_complex = pd.concat([df_combined_complex, df_AMP, df_defensive])

        # Plot the combined dataframe
        fig = px.line(df_combined_complex, x="distance", y="norm_count", animation_frame="time", color="agent", markers=True,
                      line_shape="spline")

        fig.add_vline(x=min_dist, line_width=3, line_dash="dash", line_color="orange", name="cell membrane")
        fig.add_vrect(x0=0, x1=min_dist, line_width=0, line_dash="dash", fillcolor="orange", opacity=0.2)
        fig.add_annotation(x=min_dist/2, y=df_combined_complex["norm_count"].max()/4, text="Pathogenic cell", showarrow=False)
        fig.update_layout(xaxis_title=r"$\text{Distance to cell (}\mu m\text{)}$", yaxis_title=r"$\large{\text{Concentration (agents }10^{3}/\mu m^{3})}$", legend_title="Agents", font=dict(size=15))

        with open(path+"/spatial_info_multiple_receptors.html", "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn', auto_play=False))
"""
if "bounds_uptake.csv" in os.listdir(path):
    df = pd.read_csv(path + "/bounds_uptake.csv", sep=";")
    if not(df.empty):
        df2 = pd.DataFrame()
        for k in np.unique(df["nb_bounds"]):
            d = {}
            test = []
            for i in np.unique(df["id"]):
                test.append(len(df[(df["id"] == i) & (df["nb_bounds"] == k)]))
            d["nb_bounds"] = k
            d["count"] = np.mean(test)
            d["std"] = np.std(test)
            df2 = pd.concat([df2, pd.DataFrame([d])])

        # Histogram
        fig = px.bar(df2, x="nb_bounds", y="count", error_y = "std")
        fig.update_layout(title = "AMP uptaken: nb of bounds done before", yaxis_title = "Count", xaxis_title = "nb of bounds")

        fig2 = px.box(df["nb_bounds"])
        fig2.update_layout(title = "Distribution", yaxis_title = "nb of bounds",font = dict(size = 15))
        with open(path+"/bounds_uptake.html", "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))
            f.write(fig2.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

if "bounds_distrib.csv" in os.listdir(path):
    df = pd.read_csv(path + "/bounds_distrib.csv", sep=";")
    if not df.empty:
        # Map agent values to custom labels
        df["agent_name"] = df["agent"].map({1: "AMP", 2: "Defensive"})

        df2 = pd.DataFrame()
        for j in np.unique(df["agent_name"]):
            for k in np.unique(df["nb_bounds"]):
                d = {}
                test = []
                for i in np.unique(df["id"]):
                    test.append(len(df[(df["id"] == i) & (df["agent_name"] == j) & (df["nb_bounds"] == k)]))
                d["nb_bounds"] = k
                d["agent_name"] = j
                d["count"] = np.mean(test)
                d["std"] = np.std(test)
                df2 = pd.concat([df2, pd.DataFrame([d])])

        # Histogram
        fig = px.bar(df2, x="nb_bounds", y="count", color="agent_name", opacity=0.6, barmode="overlay", error_y = "std")
        fig.update_layout(title="Nb of bounds - Histogram", yaxis_title="Count", xaxis_title="nb of bounds")

        # Box Plot
        fig2 = px.box(df, x="agent_name", y="nb_bounds", points="all", color="agent_name")
        fig2.update_layout(title="Distribution - Box Plot", yaxis_title="nb of bounds", font=dict(size=15))

        with open(path + "/bounds_distrib.html", "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))
            f.write(fig2.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

        if "bounds_uptake.csv" in os.listdir(path):
            df3 = pd.read_csv(path + "/bounds_uptake.csv", sep=";")
            if not(df3.empty):
                df4 = pd.DataFrame()
                for k in np.unique(df["nb_bounds"]):
                    d2 = {}
                    test = []
                    for i in np.unique(df3["id"]):
                        test.append(len(df3[(df3["id"] == i) & (df3["nb_bounds"] == k)]))
                    d2["nb_bounds"] = k
                    d2["count"] = np.mean(test)
                    d2["std"] = np.std(test)
                    df4 = pd.concat([df4, pd.DataFrame([d2])])
                df4["agent_name"] = "AMP_uptaken"
                df5 = pd.concat([df4, df2])
                color_map = {
                    "AMP_uptaken": "purple",
                    "AMP": "blue",
                    "Defensive": "red"
                }

                # Histogram
                fig3 = px.bar(df5, x="nb_bounds", y="count", color="agent_name", opacity=0.6, barmode="group", error_y = "std", color_discrete_map=color_map)
                fig3.update_layout(title="Nb of bounds - Histogram", yaxis_title="Count", xaxis_title="nb of bounds")
                with open(path + "/bounds_summary.html", "w") as f:
                    f.write(fig3.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))