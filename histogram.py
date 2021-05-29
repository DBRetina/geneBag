import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff

col_local_cont = list()
col_local_jac = list()
col_across_cont = list()
col_across_jac = list()

main_transcriptome_name = "ver28_GRCh37"
across_transcriptome_name = "ver38_GRCh38"


with open("ver28_GRCh37_matches.tsv", 'r') as READER:
    next(READER)
    for line in READER:
        line = line.strip()
        local, across = tuple(line.split('\t')[1:])
        # print(f"{local} | {across}")
        # continue
        if '-' in local:
            across = eval(across)
            for rec in across:
                col_across_jac.append(rec[1])
                col_across_cont.append(rec[2])

            col_local_jac.append(0)
            col_local_cont.append(0)

        elif '-' in across:
            local = eval(local)
            for rec in local:
                col_local_jac.append(rec[1])
                col_local_cont.append(rec[2])
            
            col_across_jac.append(0)
            col_across_cont.append(0)

        else:
            across = eval(across)
            local = eval(local)
            
            for rec in local:
                col_local_jac.append(rec[1])
                col_local_cont.append(rec[2])
            
            for rec in across:
                col_across_jac.append(rec[1])
                col_across_cont.append(rec[2])


df_dict_local = {'jaccard': col_local_jac, 'containment': col_local_cont}
df_dict_across = {'jaccard': col_across_jac, 'containment': col_across_cont}

df_local = pd.DataFrame(df_dict_local)
df_across = pd.DataFrame(df_dict_across)

def draw_histogram(df, title, filename):
    fig = px.histogram(df, log_y=True)
    # fig.update_layout(barmode='overlay')

    fig.update_layout(
        title_text = title,
        
        xaxis_title_text='Jac/Cont', # xaxis label
        yaxis_title_text='Count', # yaxis label
        bargap=0.2, # gap between bars of adjacent location coordinates
        bargroupgap=0.1, # gap between bars of the same location coordinates
        xaxis = dict(
            tickmode = 'linear',
        )
    )

    fig.update_traces(opacity=0.75)
    fig.write_html(f"{filename}.html")


# plot local
title_text=f"{main_transcriptome_name} local"
draw_histogram(df = df_local, title = title_text, filename=f"{main_transcriptome_name}_local")

# plot across
title_text=f"{main_transcriptome_name} across {across_transcriptome_name}" 
draw_histogram(df = df_across, title = title_text, filename=f"{main_transcriptome_name}_across_{across_transcriptome_name}")
