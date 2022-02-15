using CSV
using HTTP, JSON, DataFrames
# using PyPlot
using Plots
using RollingFunctions
using Dates
using Gumbo
using Cascadia
using Statistics
using Downloads

res = HTTP.get("https://squeezemetrics.com/monitor/static/DIX.csv")
data = CSV.read(res.body, DataFrame)

p1=plot(data[end-40:end,:price], title="SPX500", color=:green,leg=false)
p2=plot(data[end-40:end,:dix], title="DIX", color=:blue,leg=false, linestyle=:dash)
plot!([0.47], seriestype="hline", color=:black, linestyle=:dash)
p3=plot(data[end-40:end,:gex], title="GEX", color=:orange,leg=false, linestyle=:dash)
plot!([0], seriestype="hline",color=:black, linestyle=:dash)
display(plot(p1,p2,p3,layout=(3,1),legend = false))

window=50
data.gex_zscore = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[window:end,:gex_zscore] .= (data[window:end,:gex] - rollmean(data[!,:gex],window)) ./ rollstd(data[!,:gex],window)
data.dix_zscore = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[window:end,:dix_zscore] .= (data[window:end,:dix] - rollmean(data[!,:dix],window)) ./ rollstd(data[!,:dix],window)
data.price_zscore = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[window:end,:price_zscore] .= (data[window:end,:price] - rollmean(data[!,:price],window)) ./ rollstd(data[!,:price],window)
data.price_zscore.=coalesce.(data.price_zscore, NaN)
data.dix_zscore.=coalesce.(data.dix_zscore, NaN)
data.gex_zscore.=coalesce.(data.gex_zscore, NaN)

window_sma=20
data.sma_price = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[window:end,:sma_price]=rollmean(data[!,:price],window)
data.sma_dix = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[window:end,:sma_dix]=rollmean(data[!,:dix],window)
data.sma_gex = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[window:end,:sma_gex]=rollmean(data[!,:gex],window)
data.sma_dix.=coalesce.(data.sma_dix, NaN)
data.sma_gex.=coalesce.(data.sma_gex, NaN)
data.sma_price.=coalesce.(data.sma_price, NaN)

cond1=(data.dix.<0.4)
cond2=(data.dix.<data.sma_dix)
cond4=(data.gex.<data.sma_gex)
cond3=(data.price.>data.sma_price)
data.signal_sma = ones(Int, (size(data,1)))
data[cond1 .& cond2 .& cond3 .& cond4, :signal_sma] .= -1

n=400
p1 = plot(data[end-n:end, :price_zscore], title="SPX500", color=:green,leg=false)
plot!(data[end-n:end, :dix_zscore], title="DIX", color=:blue,leg=false, linestyle=:dash)
plot!([0.47], seriestype="hline", color=:black, linestyle=:dash)
p2 = plot(data[end-n:end, :price_zscore], title="SPX500", color=:green,leg=false)
plot!(data[end-n:end, :gex_zscore], title="GEX", color=:orange,leg=false, linestyle=:dash)
display(plot(p1, p2, layout = (2, 1), legend = false))

n=300
p1 = plot(data[end-n:end, :dix], title="Dix", color=:black,leg=false)
plot!([0.39], seriestype="hline", color=:red, linestyle=:dash)
plot!([0.47], seriestype="hline", color=:pink, linestyle=:dash)
p2 = plot(data[end-n:end, :price], title="Long Signals Only", color=:green,leg=false)
display(plot(p1, p2, layout = (2, 1), legend = false))

n=200
p1 = plot(data[end-n:end, :dix], title="Dix", color=:black,leg=false)
plot!(data[end-n:end, :sma_dix], color=:blue, linestyle=:dash)
plot!([0.4], seriestype="hline", color=:red, linestyle=:dash)
p2 = plot(data[end-n:end, :gex], title="Gex", color=:black,leg=false)
plot!(data[end-n:end, :sma_gex], color=:orange, linestyle=:dash)
p3 = plot(data[end-n:end, :price], title="Short Signals Only", color=:green,leg=false)
plot!(data[end-n:end, :sma_price], color=:black, linestyle=:dash)
p4=plot(data[end-n:end,:signal_sma],seriestype=:sticks, title="Signal", color=:purple,leg=false)
display(plot(p1, p2, p3,p4, layout = (4, 1), legend = false))

starti=datetime2unix(DateTime(2010,1,1))
starti=floor(Int,starti)
endi=datetime2unix(round(DateTime(today()),Dates.Day(1)))
endi=floor(Int,endi)

url=string("https://query1.finance.yahoo.com/v7/finance/download/%5EVIX?period1=",string(starti),"&period2=",string(endi),"&interval=1d&events=history")
VIX= CSV.read(Downloads.download(url),DataFrame)
rename!(VIX, Symbol("Date") => :date)
VIX[!,:VIX].=[x for x in VIX[!,:Close]]
VIX.VIX.=coalesce.(VIX.VIX,NaN)
for i = 1 : length(VIX.VIX)
    if isnothing(VIX.VIX[i])
        VIX.VIX[i] = VIX.VIX[i-1]
    end
end
VIX[!,:VIX] = [x for x in VIX.VIX]
VIX=VIX[!,[:VIX,:date]]
data=leftjoin(data,VIX,on=:date)
data.VIX_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:VIX_sma50].=rollmean(data[!,:VIX],50)
data.VIX_sma50.=coalesce.(data.VIX_sma50,NaN)
data.VIX_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:VIX_sma200].=rollmean(data[!,:VIX],200)
data.VIX_sma200.=coalesce.(data.VIX_sma200,NaN)
data.VIX_corr15 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[15:end,:VIX_corr15].=rollcor(data[!,:VIX],Vector{Union{Missing,Float64}}(data[!,:price]),15)
data.VIX_corr15.=coalesce.(data.VIX_corr15,NaN)

n=300
p1 = plot(data[end-n:end,:price], title="SPX500", color=:green,leg=false)
p2 = plot(data[end-n:end,:VIX], title="VIX", color=:black,leg=false)
p3 = plot(data[end-n:end,:VIX_corr15], title="Correlation 15D", color=:blue,leg=false)
plot!([0], seriestype="hline", color=:red, linestyle=:dash)
display(plot(p1, p2, p3, layout = (3, 1), legend = false))

url=string("https://query1.finance.yahoo.com/v7/finance/download/HYG?period1=",string(starti),"&period2=",string(endi),"&interval=1d&events=history")
HYG= CSV.read(Downloads.download(url),DataFrame)
rename!(HYG, Symbol("Date") => :date)
HYG[!,:HYG].=[x for x in HYG[!,:Close]]
HYG.HYG.=coalesce.(HYG.HYG,NaN)
for i = 1 : length(HYG.HYG)
    if isnothing(HYG.HYG[i])
        HYG.HYG[i] = HYG.HYG[i-1]
    end
end
HYG[!,:HYG] = [x for x in HYG.HYG]
HYG=HYG[!,[:HYG,:date]]
data=leftjoin(data,HYG,on=:date)
data.HYG_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:HYG_sma50].=rollmean(data[!,:HYG],50)
data.HYG_sma50.=coalesce.(data.HYG_sma50, NaN)
data.HYG_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:HYG_sma200].=rollmean(data[!,:HYG],200)
data.HYG_sma200.=coalesce.(data.HYG_sma200, NaN)

url=string("https://query1.finance.yahoo.com/v7/finance/download/LQD?period1=",string(starti),"&period2=",string(endi),"&interval=1d&events=history")
LQD= CSV.read(Downloads.download(url),DataFrame)
rename!(LQD, Symbol("Date") => :date)
LQD[!,:LQD].=[x for x in LQD[!,:Close]]
LQD.LQD.=coalesce.(LQD.LQD,NaN)
for i = 1 : length(LQD.LQD)
    if isnothing(LQD.LQD[i])
        LQD.LQD[i] = LQD.LQD[i-1]
    end
end
LQD[!,:LQD] = [x for x in LQD.LQD]
LQD=LQD[!,[:LQD,:date]]
data=leftjoin(data,LQD,on=:date)
data.LQD_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:LQD_sma50].=rollmean(data[!,:LQD],50)
data.LQD_sma50.=coalesce.(data.LQD_sma50,NaN)
data.LQD_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:LQD_sma200].=rollmean(data[!,:LQD],200)
data.LQD_sma200.=coalesce.(data.LQD_sma200,NaN)

url=string("https://query1.finance.yahoo.com/v7/finance/download/KRE?period1=",string(starti),"&period2=",string(endi),"&interval=1d&events=history")
KRE= CSV.read(Downloads.download(url),DataFrame)
rename!(KRE,Symbol("Date")=>:date)
KRE[!,:KRE].=[x for x in KRE[!,:Close]]
KRE.KRE.=coalesce.(KRE.KRE,NaN)
for i = 1 : length(KRE.KRE)
    if isnothing(KRE.KRE[i])
        KRE.KRE[i] = KRE.KRE[i-1]
    end
end
KRE[!,:KRE] = [x for x in KRE.KRE]
KRE=KRE[!,[:KRE,:date]]
data=leftjoin(data,KRE,on=:date)
data.KRE_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:KRE_sma50].=rollmean(data[!,:KRE],50)
data.KRE_sma50.=coalesce.(data.KRE_sma50,NaN)
data.KRE_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:KRE_sma200].=rollmean(data[!,:KRE],200)
data.KRE_sma200.=coalesce.(data.KRE_sma200,NaN)

n=40
p1 = plot(data[end-n:end, :price], title="SPX500", color=:green,leg=false)
p2 = plot(data[end-n:end, :HYG], title="High Yield Corporate Bond", color=:black,leg=false)
plot!(data[end-n:end, :HYG_sma50], color=:blue)
plot!(data[end-n:end, :HYG_sma200], color=:red)
p3 = plot(data[end-n:end, :LQD], title="Investment Grade Corporate Bond", color=:black,leg=false)
plot!(data[end-n:end, :LQD_sma50], color=:blue)
plot!(data[end-n:end, :LQD_sma200], color=:red)
p4 = plot(data[end-n:end, :KRE], title="Regional Banking", color=:black,leg=false)
plot!(data[end-n:end, :KRE_sma50], color=:blue)
plot!(data[end-n:end, :KRE_sma200], color=:red)
display(plot(p1, p2, p3, p4, layout = (4, 1), legend = false))

url=string("https://query1.finance.yahoo.com/v7/finance/download/CL=F?period1=",string(starti),"&period2=",string(endi),"&interval=1d&events=history")
WTIUSD = CSV.read(Downloads.download(url),DataFrame)
rename!(WTIUSD, Symbol("Date") => :date)
WTIUSD[!,:WTIUSD].=[tryparse(Float64,x) for x in WTIUSD[!,:Close]]
WTIUSD.WTIUSD.=coalesce.(WTIUSD.WTIUSD,NaN)
for i = 1 : length(WTIUSD.WTIUSD)
    if isnothing(WTIUSD.WTIUSD[i])
        WTIUSD.WTIUSD[i] = WTIUSD.WTIUSD[i-1]
    end
end
WTIUSD[!,:WTIUSD] = [x for x in WTIUSD.WTIUSD]
WTIUSD=WTIUSD[!,[:WTIUSD,:date]]
data=leftjoin(data,WTIUSD,on=:date)
data=sort(data,:date)
data.WTIUSD_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:WTIUSD_sma50].=rollmean(data[!,:WTIUSD],50)
data.WTIUSD_sma50.=coalesce.(data.WTIUSD_sma50, NaN)
data.WTIUSD_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:WTIUSD_sma200].=rollmean(data[!,:WTIUSD],200)
data.WTIUSD_sma200.=coalesce.(data.WTIUSD_sma200, NaN)

url=string("https://query1.finance.yahoo.com/v7/finance/download/HG=F?period1=",string(starti),"&period2=",string(endi),"&interval=1d&events=history")
COPPER = CSV.read(Downloads.download(url),DataFrame)
rename!(COPPER, Symbol("Date") => :date)
COPPER[!,:COPPER].=[tryparse(Float64,x) for x in COPPER[!,:Close]]
COPPER.COPPER.=coalesce.(COPPER.COPPER, NaN)
for i = 1 : length(COPPER.COPPER)
    if isnothing(COPPER.COPPER[i])
        COPPER.COPPER[i] = COPPER.COPPER[i-1]
    end
end
COPPER[!,:COPPER] = [x for x in COPPER.COPPER]
COPPER=COPPER[!,[:COPPER,:date]]
data=leftjoin(data,COPPER,on=:date)
data=sort(data,:date)
data.COPPER_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:COPPER_sma50].=rollmean(data[!,:COPPER],50)
data.COPPER_sma50.=coalesce.(data.COPPER_sma50, NaN)
data.COPPER_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:COPPER_sma200].=rollmean(data[!,:COPPER],200)
data.COPPER_sma200.=coalesce.(data.COPPER_sma200, NaN)

p1 = plot(data[end-n:end, :price], title="SPX500", color=:green,leg=false)
p2 = plot(data[end-n:end, :WTIUSD], title="WTIUSD", color=:black,leg=false)
plot!(data[end-n:end, :WTIUSD_sma50], color=:blue)
plot!(data[end-n:end, :WTIUSD_sma200], color=:red)
p3 = plot(data[end-n:end,:COPPER], title="COPPER", color=:black,leg=false)
plot!(data[end-n:end,:COPPER_sma50], color=:blue)
plot!(data[end-n:end,:COPPER_sma200], color=:red)
display(plot(p1, p2, p3, layout = (3, 1), legend = false))

resp = HTTP.get("https://www.quandl.com/api/v3/datasets/USTREASURY/YIELD.json?api_key=1y9mv-rePWeheuTZvrvx")
str = String(resp.body)
job=JSON.parse(str)
DGS30 = Array{Union{Nothing,Float64},1}(nothing,size(job["dataset"]["data"],1))
DGS2 = Array{Union{Nothing,Float64},1}(nothing,size(job["dataset"]["data"],1))
date = Array{Union{Nothing,String},1}(nothing,size(job["dataset"]["data"],1))
for i=1:size(job["dataset"]["data"],1)
    for j=1:size(job["dataset"]["data"][i])[1]
        if j==1
            date[i]=job["dataset"]["data"][i][j]
        elseif j==6
            DGS2[i]=job["dataset"]["data"][i][j]
        elseif j==12
            DGS30[i]=job["dataset"]["data"][i][j]
        end
    end
end
DGS=DataFrame(hcat(date,DGS2,DGS30),:auto)
rename!(DGS, Symbol("x1") => :date)
rename!(DGS, Symbol("x2") => :DGS2)
rename!(DGS, Symbol("x3") => :DGS30)
DGS[!,:date] = [tryparse(Date,x) for x in DGS[!,:date]]
DGS=sort(DGS,:date)
for i = 2 : length(DGS.DGS2)
    if isnothing(DGS.DGS2[i])
        DGS.DGS2[i] = DGS.DGS2[i-1]
    end
    if isnothing(DGS.DGS30[i])
        DGS.DGS30[i] = DGS.DGS30[i-1]
    end
end
replace!(DGS.DGS30,nothing=>NaN)
replace!(DGS.DGS30,missing=>NaN)
DGS[!,:DGS30] = [x for x in DGS.DGS30]
DGS[!,:DGS2] = [x for x in DGS.DGS2]
data=leftjoin(data,DGS,on=:date)
data=sort(data,:date)
data.DGS2_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:DGS2_sma50].=rollmean(data[!,:DGS2],50)
data.DGS2_sma50.=coalesce.(data.DGS2_sma50,NaN)
data.DGS2_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:DGS2_sma200].=rollmean(data[!,:DGS2],200)
data.DGS2_sma200.=coalesce.(data.DGS2_sma200,NaN)
data.DGS30_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:DGS30_sma50].=rollmean(data[!,:DGS30],50)
data.DGS30_sma50.=coalesce.(data.DGS30_sma50,NaN)
data.DGS30_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:DGS30_sma200].=rollmean(data[!,:DGS30],200)
data.DGS30_sma200.=coalesce.(data.DGS30_sma200,NaN)

p1 = plot(data[end-n:end,:price], title="SPX500", color=:green,leg=false)
p2 = plot(data[end-n:end,:DGS2], title="DGS2", color=:black,leg=false)
plot!(data[end-n:end,:DGS2_sma50], color=:blue)
plot!(data[end-n:end,:DGS2_sma200], color=:red)
p3 = plot(data[end-n:end,:DGS30], title="DGS30", color=:black,leg=false)
plot!(data[end-n:end,:DGS30_sma50], color=:blue)
plot!(data[end-n:end,:DGS30_sma200], color=:red)
display(plot(p1, p2, p3, layout=(3, 1), legend = false))

tickers=String[]
url = "https://finviz.com/screener.ashx?v=111&f=exch_nasd&o=-marketcap"
webpage = parsehtml(read(Downloads.download(url), String))
allrecords = eachmatch(Selector(".screener-link-primary"),webpage.root)
for q=1:length(allrecords)
    push!(tickers,nodeText(allrecords[q]))
end
url = "https://finviz.com/screener.ashx?v=111&f=exch_nasd&o=-marketcap&r=21"
webpage = parsehtml(read(Downloads.download(url), String))
allrecords = eachmatch(Selector(".screener-link-primary"),webpage.root)
for q=1:length(allrecords)
    push!(tickers,nodeText(allrecords[q]))
end
url = "https://finviz.com/screener.ashx?v=111&f=exch_nasd&o=-marketcap&r=41"
webpage = parsehtml(read(Downloads.download(url), String))
allrecords = eachmatch(Selector(".screener-link-primary"),webpage.root)
for q=1:length(allrecords)
    push!(tickers,nodeText(allrecords[q]))
end
url = "https://finviz.com/screener.ashx?v=111&f=exch_nasd&o=-marketcap&r=61"
webpage = parsehtml(read(Downloads.download(url), String))
allrecords = eachmatch(Selector(".screener-link-primary"),webpage.root)
for q=1:length(allrecords)
    push!(tickers,nodeText(allrecords[q]))
end
url = "https://finviz.com/screener.ashx?v=111&f=exch_nasd&o=-marketcap&r=81"
webpage = parsehtml(read(Downloads.download(url), String))
allrecords = eachmatch(Selector(".screener-link-primary"),webpage.root)
for q=1:length(allrecords)
    push!(tickers,nodeText(allrecords[q]))
end

url=string("https://query1.finance.yahoo.com/v7/finance/download/",tickers[1],"?period1=",string(starti),"&period2=",string(endi),"&interval=1d&events=history")
XXX = CSV.read(Downloads.download(url),DataFrame)
rename!(XXX, Symbol("Date") => :date)
XXX[!,Symbol(tickers[1])].=[x for x in XXX[!,:Close]]
XXX[!,Symbol(tickers[1])]=coalesce.(XXX[!,Symbol(tickers[1])],NaN)
for i = 1 : length(XXX[!,Symbol(tickers[1])])
    if isnothing(XXX[!,Symbol(tickers[1])][i])
        XXX[!,Symbol(tickers[1])][i] = XXX[!,Symbol(tickers[1])][i-1]
    end
end
XXX[!,Symbol(tickers[1])] = [x for x in XXX[!,Symbol(tickers[1])]]
XXX=XXX[!,[:date,Symbol(tickers[1])]]
tickers_df=XXX

for t in 2:length(tickers)
    local url=string("https://query1.finance.yahoo.com/v7/finance/download/",tickers[t],"?period1=",string(starti),"&period2=",string(endi),"&interval=1d&events=history")
    local XXX = CSV.read(Downloads.download(url),DataFrame)
    rename!(XXX, Symbol("Date") => :date)
    if XXX[!,:Close][1] isa Float64
        XXX[!,Symbol(tickers[t])].=[x for x in XXX[!,:Close]]
    elseif XXX[!,:Close][1] isa String
        XXX[!,Symbol(tickers[t])].=[tryparse(Float64,x) for x in XXX[!,:Close]]
    end
    XXX[!,Symbol(tickers[t])]=coalesce.(XXX[!,Symbol(tickers[t])],NaN)
    for i = 1 : length(XXX[!,Symbol(tickers[t])])
        if isnothing(XXX[!,Symbol(tickers[t])][i])
            XXX[!,Symbol(tickers[t])][i] = XXX[!,Symbol(tickers[t])][i-1]
        end
    end
    XXX[!,Symbol(tickers[t])] = [x for x in XXX[!,Symbol(tickers[t])]]
    XXX=XXX[!,[:date,Symbol(tickers[t])]]
    global tickers_df=outerjoin(tickers_df,XXX,on=:date)
end

for t in 1:length(tickers)
    tickers_df[!,Symbol(string(tickers[t],"_sma10"))] = Array{Union{Missing,Float64},1}(missing,size(tickers_df,1))
    tickers_df[10:end,Symbol(string(tickers[t],"_sma10"))].=rollmean(tickers_df[!,Symbol(tickers[t])],10)
    tickers_df[!,Symbol(string(tickers[t],"_sma10"))].=coalesce.(tickers_df[!,Symbol(string(tickers[t],"_sma10"))],NaN)
    tickers_df[!,Symbol(string(tickers[t],"_sma50"))] = Array{Union{Missing,Float64},1}(missing,size(tickers_df,1))
    tickers_df[50:end,Symbol(string(tickers[t],"_sma50"))].=rollmean(tickers_df[!,Symbol(tickers[t])],50)
    tickers_df[!,Symbol(string(tickers[t],"_sma50"))].=coalesce.(tickers_df[!,Symbol(string(tickers[t],"_sma50"))],NaN)
    tickers_df[!,Symbol(string(tickers[t],"_sma200"))] = Array{Union{Missing,Float64},1}(missing,size(tickers_df,1))
    tickers_df[200:end,Symbol(string(tickers[t],"_sma200"))].=rollmean(tickers_df[!,Symbol(tickers[t])],200)
    tickers_df[!,Symbol(string(tickers[t],"_sma200"))].=coalesce.(tickers_df[!,Symbol(string(tickers[t],"_sma200"))],NaN)
    tickers_df[!,Symbol(string(tickers[t],"_NDXA10R"))].=tickers_df[!,Symbol(string(tickers[t],"_sma10"))].<tickers_df[!,Symbol(tickers[t])]
    tickers_df[!,Symbol(string(tickers[t],"_NDXA50R"))].=tickers_df[!,Symbol(string(tickers[t],"_sma50"))].<tickers_df[!,Symbol(tickers[t])]
    tickers_df[!,Symbol(string(tickers[t],"_NDXA200R"))].=tickers_df[!,Symbol(string(tickers[t],"_sma200"))].<tickers_df[!,Symbol(tickers[t])]
end

NDXA10R=copy(tickers_df)
NDXA50R=copy(tickers_df)
NDXA200R=copy(tickers_df)
for i in axes(tickers_df,2)
    if i>1
        if !occursin("NDXA10R",string(names(tickers_df)[i]))
            select!(NDXA10R,Not(names(tickers_df)[i]))
        end
        if !occursin("NDXA50R",string(names(tickers_df)[i]))
            select!(NDXA50R,Not(names(tickers_df)[i]))
        end
        if !occursin("NDXA200R",string(names(tickers_df)[i]))
            select!(NDXA200R,Not(names(tickers_df)[i]))
        end
    end
end
NDXA10R[!,2:end].=coalesce.(NDXA10R[!,2:end],NaN)
NDXA10R[!,2:end].=Float64.(NDXA10R[!,2:end])
NDXA50R[!,2:end].=coalesce.(NDXA50R[!,2:end],NaN)
NDXA50R[!,2:end].=Float64.(NDXA50R[!,2:end])
NDXA200R[!,2:end].=coalesce.(NDXA200R[!,2:end],NaN)
NDXA200R[!,2:end].=Float64.(NDXA200R[!,2:end])

NDXA10R[!,:NDXA10R] = fill(NaN64, size(NDXA10R)[1])
NDXA50R[!,:NDXA50R] = fill(NaN64, size(NDXA50R)[1])
NDXA200R[!,:NDXA200R] = fill(NaN64, size(NDXA200R)[1])
for i=1:size(NDXA10R)[1]
    v=Vector{Union{Missing,Float64}}(NDXA10R[i,2:end])
    NDXA10R[i,:NDXA10R]=mean(filter(!isnan,v))
    v=Vector{Union{Missing,Float64}}(NDXA50R[i,2:end])
    NDXA50R[i,:NDXA50R]=mean(filter(!isnan,v))
    v=Vector{Union{Missing,Float64}}(NDXA200R[i,2:end])
    NDXA200R[i,:NDXA200R]=mean(filter(!isnan,v))
end

NDXA10R = DataFrame(NDXA10R = NDXA10R[!,:NDXA10R])
NDXA10R.date=tickers_df[!,:date]
NDXA50R = DataFrame(NDXA50R = NDXA50R[!,:NDXA50R])
NDXA50R.date=tickers_df[!,:date]
NDXA200R = DataFrame(NDXA200R = NDXA200R[!,:NDXA200R])
NDXA200R.date=tickers_df[!,:date]

data=innerjoin(data,NDXA10R,on=:date)
data.NDXA10R_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:NDXA10R_sma50].=rollmean(data[!,:NDXA10R],50)
data.NDXA10R_sma50.=coalesce.(data.NDXA10R_sma50,NaN)
data.NDXA10R_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:NDXA10R_sma200].=rollmean(data[!,:NDXA10R],200)
data.NDXA10R_sma200.=coalesce.(data.NDXA10R_sma200,NaN)

data=innerjoin(data,NDXA50R,on=:date)
data.NDXA50R_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:NDXA50R_sma50].=rollmean(data[!,:NDXA50R],50)
data.NDXA50R_sma50.=coalesce.(data.NDXA50R_sma50,NaN)
data.NDXA50R_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:NDXA50R_sma200].=rollmean(data[!,:NDXA50R],200)
data.NDXA50R_sma200.=coalesce.(data.NDXA50R_sma200,NaN)

data=innerjoin(data,NDXA200R,on=:date)
data.NDXA200R_sma50 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[50:end,:NDXA200R_sma50].=rollmean(data[!,:NDXA200R],50)
data.NDXA200R_sma50.=coalesce.(data.NDXA200R_sma50,NaN)
data.NDXA200R_sma200 = Array{Union{Missing,Float64},1}(missing,size(data,1))
data[200:end,:NDXA200R_sma200].=rollmean(data[!,:NDXA200R],200)
data.NDXA200R_sma200.=coalesce.(data.NDXA200R_sma200,NaN)

p1 = plot(data[end-n:end, :price], title="SPX500", color=:green,leg=false)
p2 = plot(data[end-n:end, :NDXA10R], title="NDXA10R", color=:black,leg=false)
plot!(data[end-n:end, :NDXA10R_sma50], color=:blue)
plot!(data[end-n:end, :NDXA10R_sma200], color=:red)
p3 = plot(data[end-n:end, :NDXA50R], title="NDXA50R", color=:black,leg=false)
plot!(data[end-n:end, :NDXA50R_sma50], color=:blue)
plot!(data[end-n:end, :NDXA50R_sma200], color=:red)
p4 = plot(data[end-n:end, :NDXA200R], title="NDXA200R", color=:black,leg=false)
plot!(data[end-n:end, :NDXA200R_sma50], color=:blue)
plot!(data[end-n:end, :NDXA200R_sma200], color=:red)
display(plot(p1, p2, p3, p4, layout=(4, 1), legend = false))

resp = HTTP.get("https://www.quandl.com/api/v3/datasets/AAII/AAII_SENTIMENT.json?api_key=1y9mv-rePWeheuTZvrvx")
str = String(resp.body)
job=JSON.parse(str)
AAII=DataFrame(job["dataset"]["data"],:auto)
AAII=DataFrame([[names(AAII)]; collect.(eachrow(AAII))], [:column; Symbol.(axes(AAII, 1))])
select!(AAII,:2,:8)
rename!(AAII, Symbol("1")=>:date)
rename!(AAII, Symbol("7")=>:BuBeSpread)
AAII.date = Date.(AAII.date, "yyyy-mm-dd")
AAII.BuBeSpread.=coalesce.(AAII.BuBeSpread,NaN)
AAII=sort(AAII,:date)
data=leftjoin(data, AAII, on = :date)
for i=2:size(data,1)
    if ismissing(data[i,:BuBeSpread])
        data[i,:BuBeSpread]=data[i-1,:BuBeSpread]
    end
end

p1 = plot(data[end-n:end,:price], title="SPX500", color=:green,leg=false)
p2 = plot(data[end-n:end,:BuBeSpread], title="BuBeSpread", color=:brown,leg=false)
plot!([-0.2], seriestype="hline", color=:red, linestyle=:dash)
plot!([0], seriestype="hline", color=:black, linestyle=:dash)
display(plot(p1, p2, layout=(2, 1), legend = false))

resp = HTTP.get("https://www.quandl.com/api/v3/datasets/FRED/T10YIE.json?api_key=1y9mv-rePWeheuTZvrvx")
str = String(resp.body)
job=JSON.parse(str)
T10YIE=DataFrame(job["dataset"]["data"],:auto)
T10YIE=DataFrame([[names(T10YIE)]; collect.(eachrow(T10YIE))], [:column; Symbol.(axes(T10YIE, 1))])
select!(T10YIE, :2, :3)
rename!(T10YIE, Symbol("1") => :date)
rename!(T10YIE, Symbol("2") => :T10YIE)
T10YIE.date = Date.(T10YIE.date, "yyyy-mm-dd")
T10YIE.T10YIE.=coalesce.(T10YIE.T10YIE,NaN)
T10YIE=sort(T10YIE,:date)
data=leftjoin(data,T10YIE,on=:date)
data=sort(data,:date)
for i=2:size(data,1)
    if ismissing(data[i,:T10YIE])
        data[i,:T10YIE]=data[i-1,:T10YIE]
    end
end

resp = HTTP.get("https://www.quandl.com/api/v3/datasets/FRED/T5YIFR.json?api_key=1y9mv-rePWeheuTZvrvx")
str = String(resp.body)
job=JSON.parse(str)
T5YIFR=DataFrame(job["dataset"]["data"],:auto)
T5YIFR=DataFrame([[names(T5YIFR)]; collect.(eachrow(T5YIFR))], [:column; Symbol.(axes(T5YIFR, 1))])
select!(T5YIFR, :2, :3)
rename!(T5YIFR, Symbol("1") => :date)
rename!(T5YIFR, Symbol("2") => :T5YIFR)
T5YIFR.date = Date.(T5YIFR.date, "yyyy-mm-dd")
T5YIFR.T5YIFR.=coalesce.(T5YIFR.T5YIFR,NaN)
T5YIFR=sort(T5YIFR,:date)
data=leftjoin(data,T5YIFR,on=:date)
data=sort(data,:date)
for i=2:size(data,1)
    if ismissing(data[i,:T5YIFR])
        data[i,:T5YIFR]=data[i-1,:T5YIFR]
    end
end

p1 = plot(data[end-n:end,:price], title="SPX500", color=:green,leg=false)
p2 = plot(data[end-n:end,:T10YIE], title="T10YIE", color=:red,leg=false)
p3 = plot(data[end-n:end,:T5YIFR], title="T5YIFR", color=:brown,leg=false)
display(plot(p1,p2,p3, layout=(3, 1), legend = false))

data=sort(data,:date,rev=true)
