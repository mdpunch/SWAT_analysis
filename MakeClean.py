# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python (ctapipe_0.24)
#     language: python
#     name: ctapipe_0.24
# ---

# %%
MC = "MC_8MSTs_Proton_with_mono_lapalma.txt"
MCo =  "MC_8MSTs_Proton_with_mono_removed_lapalma.txt"
MC = "MC_8MSTs_Gamma_spec_2.7_with_mono_lapalma.txt"
MCo = "MC_8MSTs_Gamma_spec_2.7_with_mono_removed_lapalma.txt"
MC_file = open(MC,"r")
MC_out = open(MCo,"wt")



# %%
# Go through lines and check for >= 2 telescopes
for evt in MC_file:
    if len(evt)-1-evt.count('.') > 1:
        MC_out.write(evt)


# %%
