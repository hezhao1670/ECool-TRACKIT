using Plots,DelimitedFiles
using LaTeXStrings
#pyplot()

dir=@__DIR__
file1 = readdlm(dir*"//csmon.out")

speedfactor = 6
file = file1

kturns_column = 1
emittx_column = 2
emitty_column = 3
sigtime_column = 4
sigdpp_column = 5

xavg_column = 6
yavg_column = 7
dppavg_column =8

np_integrate = 9
np_losst = 10
np_lossinj = 11

tunex=3.62
tuney=2.63
circ=161
gamma = 1.0042
beta = sqrt(1-1/gamma^2)
f0 = beta*3e8/circ
betax = circ/(2*pi*tunex)
betay = circ/(2*pi*tuney)

turns = getindex(file,:,kturns_column)./f0*speedfactor

y1 = getindex(file,:,emittx_column)*10^6
y2 = getindex(file,:,emitty_column)*10^6
y3 = getindex(file,:,sigdpp_column)

y4 = getindex(file,:,xavg_column)
y5 = getindex(file,:,yavg_column)
y6 = getindex(file,:,dppavg_column)

y7 = getindex(file,:,np_integrate)
y8 = getindex(file,:,np_losst)

pic1=plot(turns,y2,markersize=0.2,label=L"Ver.",xlabel="Time (s)",lw=2,ls=:solid)
pic2=plot!(turns,y1,markersize=0.3,label=L"Hori.",frame=:box,
    lw=2,ylabel="Emittance (um)",ls=:solid)
pic3=plot(turns,y3,markersize=0.3,label=L"dp/p",frame=:box,
    ylabel="dp/p",dpi=250,lw=2,ls=:solid)
plot(pic2,pic3,layout=(2,1))


pic4=plot(turns,y4,markersize=0.2,label=L"x-avg..",xlabel="Time (s)",lw=2,ls=:solid)
pic5=plot!(turns,y5,markersize=0.3,label=L"y-avg.",frame=:box,
    lw=2,ylabel="Position offset(m)",ls=:solid)
pic6=plot(turns,y6,markersize=0.3,label=L"dp/p-avg",frame=:box,
    ylabel="dp/p offset",dpi=250,lw=2,ls=:solid,ylims=(-0.0002,0.0024))
plot(pic5,pic6,layout=(2,1))


pic1=plot(turns,y7,markersize=0.2,label=L"np.",xlabel="Time (s)",lw=2,frame=:box,ylims=(0,4e9))

plot(pic2,pic3,pic6,layout=(3,1),size=(500,500))
#plot(pic1,pic3,layout=(2,1))
#savefig(dir*"/evo")




function getdistri(p_data,rms_x,k_span,bins,distri)
# give the beam distribution in x axis
    width = k_span*rms_x
    bins_width = width/bins
    for i = 1:length(distri[:,1])
        distri[i,2] = 0.0
    end
    for i = 1:bins
        posi = -width/2 + (i-1)*bins_width
        distri[i,1] = posi
    end
    for i = 1:length(p_data[:,1])
        xx = p_data[i,1]
        x_posi = xx/(width/2.) * bins/2. + bins/2.
        x_posi = round(x_posi)
        x_posi = convert(Int64,x_posi)
        if x_posi<1
            x_posi = 1
        elseif x_posi>bins
            x_posi = bins
        end
        distri[x_posi,2] = distri[x_posi,2] + 1
    end
    #distri[:,2] = distri[:,2]./length(p_data[:,1])
end


##=====================================================
dir=@__DIR__
file31 = readdlm(dir*"//ini.full")
file32 = readdlm(dir*"//tran.full")
n1 = 1
n2 = 6

xxscale = 1
xi = getindex(file31,:,n1)#.*f0*2*pi/2
xpi = getindex(file31,:,n2)
xf = getindex(file32,:,n1)#.*f0*2*pi/2
xpf = getindex(file32,:,n2)


cur_colors = get_color_palette(:auto,3)

bins = 180
hist_x = zeros(Float64,bins,2)
getdistri(xi*xxscale,0.00435,8,bins,hist_x)
plot(hist_x[:,1],hist_x[:,2]./380000 .- 1.5e-3,label="",line=:solid,fill=(0.0015,0.5,:gray),lc=:RED)

hist_x = zeros(Float64,bins,2)
getdistri(xf*xxscale,0.00435,8,bins,hist_x)
plot!(hist_x[:,1],hist_x[:,2]./380000 .- 1.5e-3,label="",line=:solid,fill=(0.0015,0.5,:gray),lc=:purple)

scatter(xi*xxscale,xpi,markersize=1.15,markerstrokewidth=0,xlabel="Position (μs)",frame=:box, ylabel="dp/p",
    ylims=(-0.0016,0.0016),label="t=0.0 s",mcolor=cur_colors[2])
scatter!(xf*xxscale,xpf,markersize=1.15,markerstrokewidth=0,frame=:box,label="t=end",title="",
    dpi=100,mcolor=cur_colors[1],xlabel="Position (μs)", ylabel="dp/p")

savefig(dir*"//x-dpp")




vline!([1.0/f0/2.0*xxscale,-1.0/f0/2.0*xxscale],line=(:dash,:black),label=:none)
vline!([elength/2.0*xxscale,-elength/2.0*xxscale],line=(:dash,:black),label=:none,dpi=200,
    xlims=(-0.55,0.55),ylims=(-0.0016,0.0016),size=(420,320),title="(a)")



bins = 280
hist_x = zeros(Float64,bins,2)
getdistri(xpi,5e-4,8,bins,hist_x)
plot!(hist_x[:,2]./3000 .- 1.0/f0/2.0*xxscale,hist_x[:,1],label="",line=:shortdash,fill=(0.5,0.5,:gray),lc=:purple)

hist_x = zeros(Float64,bins,2)
getdistri(xpf,5e-4,8,bins,hist_x)
plot!(hist_x[:,2]./3000 .- 1.0/f0/2.0*xxscale,hist_x[:,1],label="",line=:shortdash,fill=(0.002,0.5,:gray),lc=:purple)



stephist(xi./1e-6,norm=false,label="t=0 s",bins=100,ylabel="Density (arb. unit)",lw=2)
fig2 = stephist!(xf./1e-6,norm=false,label="t=1.01 s",xlabel="Position (μs)",bins=100,frame=:box,dpi=200,lw=2)
#plot(fig1,fig2)




##==== IBS =================================================
dir=@__DIR__
file00 = readdlm(dir*"\\800U92-vkick240.5-1cm-e-1store//fort.66")
n1 = 1
taux = 5
tauy = 6
tauz = 7


turns = getindex(file00,:,n1)./f0*speedfactor
taux = getindex(file00,:,5)
tauy = getindex(file00,:,6)
taup = getindex(file00,:,7)
plot(turns,[taux tauy taup],markersize=2.25,markerstrokewidth=0)



## mountain plot
dir=@__DIR__
file41 = readdlm(dir*"\\800U92-vkick240.5-gap1s21//fort.62")

span = 101
scale = 10000
length(file41[:,1])/span

data = file41[1:span,:]
data2 = data[3:span-1,:]
t = getindex(data2,:,1).*1e6
d = getindex(data2,:,2)
turns = getindex(data2,:,5)
d = d*5000
mountain = plot(t,d*scale.*1.06e-6,label=:none,frame=:box,xlabel="Time (μs)",lc=:black)

cur_colors = get_color_palette(:auto,3)

for k=1:3:457
    data = file41[1+k*span:(k+1)*span,:]
    data2 = data[3:span-1,:]
    t = getindex(data2,:,1)*1e6
    d = getindex(data2,:,2)
    turns = getindex(data2,:,5)

    d = d +turns/scale
    mountain = plot!(t,d*scale.*(speedfactor/f0),label=:none,frame=:box,xlabel="Time (μs)",ylabel="Time (s)",lc=:black)
end
mountain




## heatmap plot
dir=@__DIR__
file41 = readdlm(dir*"\\800U92-vkick240.5-gap1s21//fort.62")

file_out = [ "time"  "t-dis"   "data"]

span = 101
scale = Int(size(file41)[1]/span)-1

for k=0:scale
    data = file41[1+k*span:(k+1)*span,:]
    data2 = data[3:span-1,:]
    t_dis = getindex(data2,:,1)*1e6
    d = getindex(data2,:,2)
    turns = getindex(data2,:,5).*(speedfactor/f0)

    data_new = hcat(turns,t_dis,d)
    global file_out = vcat(file_out, data_new)
end


open(dir*"/mountain_dis.out","w") do io
    writedlm(io,file_out)
end





## heatmap plot-2
dir=@__DIR__
file41 = readdlm(dir*"\\800U92-vkick240.5-gap1s21//fort.62")

span = 101
scale = Int(size(file41)[1]/span)-1

data = file41[1:span,:]
data2 = data[3:span-1,:]
t_dis = getindex(data2,:,1)*1e6
d = getindex(data2,:,2)
time_p = getindex(data2,:,5).*(speedfactor/f0)

t_dis_out = vcat("posi",t_dis)
data_out = vcat(time_p[1],d./(floor(time_p[1])+1))
data_new = hcat(t_dis_out,data_out)


for k=1:scale
    data = file41[1+k*span:(k+1)*span,:]
    data2 = data[3:span-1,:]
    t_dis = getindex(data2,:,1)*1e6
    d = getindex(data2,:,2)
    #d = getindex(data2,:,2)./((floor(getindex(data2,:,5)[1].*(speedfactor/f0))+1.0)*0.8)
    time_p = getindex(data2,:,5).*(speedfactor/f0)

    data_out = vcat(time_p[1],d)
    data_new = hcat(data_new,data_out)
end

open(dir*"/map_dis.out","w") do io
    writedlm(io,data_new)
end

start_i = 2
end_i = scale
data = data_new[2:99,start_i:end_i]

fig_map = heatmap(data_new[1,start_i:end_i],data_new[2:99,1], data,c=:jet1,frame=:box,
    xlabel="Time (s)", ylabel="Position (μs)",colorbar=:none,colorbar_scale=:none)

#contour(1:size(data,2),1:size(data,1), data,fill=true)


pic3=plot(turns,y3.*0.1e4,markersize=0.3,label="10³δₚ",frame=:box,ylabel="dp/p",
    dpi=250,lw=2,ls=:solid,ylims=(0.13,0.46),xlims=(0,data_new[1,end_i]),grid=:true)

plot(pic3,fig_map,layout=grid(2,1,heights=[0.47,0.53]),dpi=200,size=(500,400))
#savefig(dir*"/colormap_zoomin")
