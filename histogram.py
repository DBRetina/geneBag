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
        if '-' in local:
            across = eval(across)
            across_jac, across_cont = across[0][1], across[0][2]
            col_across_cont.append(across_cont)
            col_across_jac.append(across_jac)
        elif '-' in across:
            local = eval(local)
            local_jac, local_cont = local[0][1], local[0][2]
            col_local_cont.append(local_cont)
            col_local_jac.append(local_jac)
        else:
            across = eval(across)
            local = eval(local)
            across_cont, across_jac = across[0][1], across[0][2]
            local_jac, local_cont = local[0][1], local[0][2]
            col_across_cont.append(across_cont)
            col_across_jac.append(across_jac)
            col_local_cont.append(local_cont)
            col_local_jac.append(local_jac)


df_dict_local = {'local_jac': col_local_jac, 'local_cont': col_local_cont}
df_dict_across = {'across_jac': col_across_jac, 'across_cont': col_across_cont}

df_local = pd.DataFrame(df_dict_local)
fig = px.histogram(df_local, log_y=True)
fig.update_layout(barmode='overlay')

fig.update_layout(
    # title_text=f"{main_transcriptome_name} across {across_transcriptome_name}", # title of plot
    title_text=f"{main_transcriptome_name} local", # title of plot
    
    xaxis_title_text='Jac/Cont', # xaxis label
    yaxis_title_text='Count', # yaxis label
    bargap=0.2, # gap between bars of adjacent location coordinates
    bargroupgap=0.1, # gap between bars of the same location coordinates
    xaxis = dict(
        tickmode = 'linear',
    )
)

fig.update_traces(opacity=0.75)
fig.write_html(f"{main_transcriptome_name}_local.html")

df_across = pd.DataFrame(df_dict_across)
fig = px.histogram(df_across, log_y=True)
fig.update_layout(barmode='overlay')

fig.update_layout(
    title_text=f"{main_transcriptome_name} across {across_transcriptome_name}", # title of plot
    xaxis_title_text='Jac/Cont', # xaxis label
    yaxis_title_text='Count', # yaxis label
    bargap=0.2, # gap between bars of adjacent location coordinates
    bargroupgap=0.1, # gap between bars of the same location coordinates
    xaxis = dict(
        tickmode = 'linear',
    )
)

fig.update_traces(opacity=0.75)
fig.write_html(f"{main_transcriptome_name}_across_{across_transcriptome_name}.html")
