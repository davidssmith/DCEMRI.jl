using PyPlot
using DCEMRI

rt = readcsv("runtimes.csv")
nt = float(rt[:,2])
ft = float(1./rt[:,3])

f(x::Vector{Float64}, a::Vector{Float64}) = a[1]*x.^a[2]
p0 = [1e-5, 2.0]
params, resid = nlsfit(f, ft.'.', [1], nt, p0)
println("fit params")
println(params)

figure(1, figsize=(6,5))
clf()
loglog(nt, ft, "ko")
x = logspace(0,log(maximum(nt)), 100)
loglog(x, f(x, params[:]), "k:", alpha=0.5)
xlabel("number of time points N in tissue curve")
ylabel("time to fit one voxel (s)")
ylim(1e-4,1)
xlim(10,1e4)
legend(["measured", @sprintf "%.1eN^%.1f" params[1] params[2]],
	loc="upper left", frameon=false, numpoints=1, fontsize=10)
savefig("rta.pdf")


