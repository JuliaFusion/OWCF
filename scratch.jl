### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ b910b53a-a83c-11f1-1c4d-290ec76846d3
using Pkg; Pkg.activate("."); using JLD2; using Interpolations; using Plots

# ╔═╡ 29c89906-87d4-463f-924b-7df6bd2ee15c
begin
	folderpath_OWCF = "REPLACE-THIS-TEXT-WITH-THE-PATH-TO-THE-OWCF-FOLDER-ON-YOUR-COMPUTER"
myfile = jldopen(folderpath_OWCF*"apps/example_data/ps2WF_results_JET_99971L72_at8,9s_AB_T-D--n-4He_126x101x104.jld2", false, false, false, IOStream)
	q = Dict()
	for key in keys(myfile)
		q[key] = myfile[key]
	end
	close(myfile)
end

# ╔═╡ e8cf668b-0173-4e77-9d19-27635cc44a2a


# ╔═╡ Cell order:
# ╠═b910b53a-a83c-11f1-1c4d-290ec76846d3
# ╠═29c89906-87d4-463f-924b-7df6bd2ee15c
# ╠═e8cf668b-0173-4e77-9d19-27635cc44a2a
