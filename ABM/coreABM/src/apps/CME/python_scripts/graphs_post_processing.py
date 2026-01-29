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

        num_compartments = 25
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
        df_test = df.groupby(["id", "time", "agent", "compartment"], as_index=False)["count"].sum()
        df_test["count"] = df_test["count"]/1000
        df_test["norm_count"] = df_test.apply(lambda row: row["count"] / (compartment_volume[row["compartment"]]), axis=1)
        df_res = df_test.groupby(["time", "agent", "compartment"], as_index=False)["norm_count"].mean()
        if (len(np.unique(df["id"])) > 1):
            std = df_test.groupby(["time", "agent", "compartment"], as_index=False)["norm_count"].std()["norm_count"].tolist()
        else:
            std = np.zeros(df_res.shape[0])
        df_res["std"] = std
        df_res["std"] = df_res["std"].fillna(0)

        df_res["distance"] = df_res["compartment"]*step + min_dist
        df_res = df_res.where(df_res["agent"] != "Cell").dropna()
        fig = px.line(df_res, x = "distance", y = "norm_count",
              error_y = "std", animation_frame="time", color="agent", markers = True, 
              line_shape="spline", color_discrete_map={
            "AMP": "blue",
            "Complex": "green",
            "Defensive": "red"
        })

        fig.add_vline(x=min_dist, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")
        fig.add_vrect(x0=0, x1=min_dist, line_width=0, line_dash = "dash", fillcolor="orange", opacity=0.2)
        fig.add_annotation(x = min_dist/2, y = df_test["norm_count"].max()/4, text = "Pathogenic cell", showarrow=False)
        fig.update_layout(xaxis_title = r"$\text{Distance to cell (}\mu m\text{)}$", yaxis_title = r"$\large{\text{Concentration (}10^{3} \text{agents}/\mu m^{3})}$", legend_title = "Agents", font = dict(size = 15))

        with open(path+"/post_processing_errorbars.html", "w") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn',  auto_play=False))

        fig = px.line(df_res, x = "distance", y = "norm_count",
              animation_frame="time", color="agent", markers = True, 
              line_shape="spline", color_discrete_map={
            "AMP": "blue",
            "Complex": "green",
            "Defensive": "red"
        })

        fig.add_vline(x=min_dist, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")
        fig.add_vrect(x0=0, x1=min_dist, line_width=0, line_dash = "dash", fillcolor="orange", opacity=0.2)
        fig.add_annotation(x = min_dist/2, y = df_test["norm_count"].max()/4, text = "Pathogenic cell", showarrow=False)
        fig.update_layout(xaxis_title = r"$\text{Distance to cell (}\mu m\text{)}$", yaxis_title = r"$\large{\text{Concentration (}10^{3} \text{agents}/\mu m^{3})}$", legend_title = "Agents", font = dict(size = 15))

        for nb_frame in [0, 10, 25, 50, 100, 200, 500]:
            if nb_frame < len(fig.frames):
                frame_data = fig.frames[nb_frame].data

                # Create a new Figure object with only the desired frame
                fig_single_frame = go.Figure(data=frame_data)
                fig_single_frame.add_vline(x=min_dist, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")
                fig_single_frame.add_vrect(x0=0, x1=min_dist, line_width=0, line_dash = "dash", fillcolor="orange", opacity=0.2)
                fig_single_frame.add_annotation(x = min_dist/2, y = df_test["norm_count"].max()/4, text = "Pathogenic cell", showarrow=False)
                fig_single_frame.update_layout(xaxis_title = r"$\text{Distance to cell (}\mu m\text{)}$", yaxis_title = r"$\large{\text{Concentration (}10^{3} \text{agents}/\mu m^{3})}$", legend_title = "Agents", font = dict(size = 15))

                df_res_tmp = df_res.where(df_res["time"] == df_res["time"].unique()[nb_frame]).dropna()
                count = 0

                colors = ['rgba(31, 119, 180, 0.3)', 'rgba(44, 160, 44, 0.3)', 'rgba(214, 39, 40, 0.3)']
                for i in np.unique(df_res_tmp["agent"]):
                        df_tmp = df_res_tmp.where(df_res_tmp["agent"] == i).dropna()
                        
                        fig_single_frame.add_trace(go.Scatter(x=df_tmp['distance'], y=df_tmp["norm_count"] + df_tmp["std"],
                                                 line=dict(shape = "spline", width=0, color = 'rgba(0, 0, 0, 0)'),
                                                 showlegend=False))
                        fig_single_frame.add_trace(go.Scatter(x=df_tmp['distance'], y=df_tmp["norm_count"] - df_tmp["std"],
                                                 line=dict(shape = "spline", width=0, color = 'rgba(0, 0, 0, 0)'),
                                                 fill='tonexty',
                                                 fillcolor=colors[count],
                                                 showlegend=False))
                        count += 1

                time = str(df_res["time"].unique()[nb_frame])
                fig_single_frame.add_annotation(xref="paper", yref="paper", x = 1.093, y = 0.5, text = "t = " + time +" s", showarrow=False, font=dict(size=18))
                fig_single_frame.update_yaxes(range = [-0.01, 0.5])

                with open(path+"/post_processing_errorbars.html", "a") as f:
                    f.write(fig_single_frame.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn',  auto_play=False))