import streamlit as st
from gifs.get_gifs import create_gif, get_relative_plot

PROTEIN_DICTIONARY = {
    'ORF1AB': (266, 21555),
    'S': (21563, 25384),
    'ORF3a': (25393, 26220),
    'E': (26245, 26472),
    'M': (26523, 27191),
    'ORF6': (27202, 27387),
    'ORF7a': (27394, 27759),
    'ORF8': (27894, 28259),
    'N': (28274, 29533),
    'ORF10': (29558, 29674),
    'complete_genome': (0, 30000),
}

gene = st.sidebar.selectbox(
    "Select the software:",
    options=list(PROTEIN_DICTIONARY.keys()),
)
# if apobec == "snippy":
#     df = pd.read_csv("./software/snippy/monkeypox_data.csv")
# else:
#     df = pd.read_csv("./software/ivar/monkeypox_data.csv")
#
# st.sidebar.header("Please filter Here:")
# dataset = st.sidebar.multiselect(
#     "Select the Datasets:", options=df["dataset"].unique(), default=df["dataset"].unique()
# )
# type = st.sidebar.multiselect(
#     "Select the Type:", options=df["type"].unique(), default=df["type"].unique()
# )
#
# annotation = st.sidebar.multiselect(
#     "Select the Annotation:",
#     options=df["mutation_anotation"].unique(),
#     default=df["mutation_anotation"].unique()
# )
# apobec = st.sidebar.multiselect(
#     "Select the Apobec:",
#     options=df["apobec"].unique(),
#     default=df["apobec"].unique()
# )
bins = st.sidebar.slider('Number of bins', 0.0, 1.0, 0.1)
# default_value = 0.0
max_freq_value = st.sidebar.slider('Max Freq', 0.0, 1.0,  value=1.0)
duration = st.sidebar.slider('duration (s)', 0.0, 10.0, value=1.0)
# limit = st.sidebar.slider('Limit', 1.0, 0.0)


tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(
    ["Analyse 40k Sars", "Snps prediction", "Artificial Reads", "Artificial Tests", "Real Comparison","Specific Cases"])
with tab1:
    if gene:
        get_relative_plot(PROTEIN_DICTIONARY[gene], gene)
        st.image(f"./gifs/images/relative_frequency_{gene}.png")
        create_gif(gene, bins, PROTEIN_DICTIONARY[gene], duration,
                   (0, max_freq_value))
        st.image(
            f"./gifs/images/{gene}_{bins}_{duration}.gif", use_column_width=True)
with tab2:
    st.image("./prediction/prediction.png")

with tab3:
    st.image("./artificial_data/Artifical_Data_50bins_100000_samples.png")
    st.image("./artificial_data/Pos2022_2023.png")
with tab4:

    st.image("./artificial_data/acc_snps.png")
    st.image("./artificial_data/MM_snps.png")
    st.image("./artificial_data/Ns_snps.png")
    st.image("./artificial_data/acc_identity.png")
    st.image("./artificial_data/MM_identity.png")
    st.image("./artificial_data/Ns_identity.png")
with tab5:
    st.image("./real_data/boxplot_alignment_problem.png")
    st.image("./real_data/boxplot_matches.png")
    st.image("./real_data/boxplot_gap_ivar.png")
    st.image("./real_data/boxplot_snippy_conservative.png")
    st.image("./real_data/boxplot_gap_snippy.png")
    st.image("./real_data/boxplot_snippy_ivar_mismatch.png")
    st.image("./real_data/boxplot_ivar_conservative.png")
    st.image("./real_data/boxplot_uncovered_case.png")
with tab6:
    st.image("./real_data/s_case_1", use_column_width=True)
    st.image("./real_data/s_case_2", use_column_width=True)
    st.image("./real_data/s_case_3", use_column_width=True)
    st.image("./real_data/s_case_4", use_column_width=True)
