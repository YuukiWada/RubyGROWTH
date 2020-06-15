#!/usr/bin/env ruby
# coding: utf-8
require "RubyROOT"
require "RubyFits"
require "time"
include Root
include Fits
include Math

class Growth
  def initialize(input_file)
    @fits_ins=Array.new
    if input_file.class="Array" then
      input_file.each do |file|
        @fits_ins << Fits::FitsFile.new(file)
      end
    elsif input_file.class="String"
      @fits_ins << Fits::FitsFile.new(input_file)
    else
      "Error: Input file names should be given by Array or String."
      exit 1
    end

    # warning
    # 連続したファイルでないと警告を出すようにする
    
  end

  def unixtime(string, zone="Asia/Tokyo")
    # https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
    ENV["TZ"]=zone
    time=Time.parse(string).to_i.to_f
    return time
  end

  def define_bins(min, max, num)
    log_min=log(min, 10)
    log_max=log(max, 10)
    width=(log_max-log_min)/num.to_f
    bins=Root::DoubleArray.new(num+1)
    (num+1).times do |nbin|
      bins[nbin]=10**(log_min+width*nbin.to_f)
    end
    bins[0]=min
    bins[num]=max
    return bins
  end

  def ql_spec(adc=0, min=0.05, max=50.0, num=200, log=true, name="ql_spec")
    if @fits_ins.length!=1 then
      puts "Error: Quick look mode can be used with a single fits file."
      exit 1
    end
    
    if log==true then
      bins=define_bins(min, max, num)
    else
      bins=Root::DoubleArray.new(num+1)
      width=(max-min)/num.to_f
      for i in 0..num
        bins[i]=min+width*i.to_f
      end
      bins[0]=min
      bins[num]=max
    end
    
    hist=Root::TH1D.create(name, name, num, bins)
    fits=@fits_ins[0]
    event_hdu=fits[0]["EVENTS"]
    event_num=event_hdu.getNRows()-1
    adc_index=event_hdu["boardIndexAndChannel"]
    energy_raw=event_hdu["energy"]
    energy_width_header="BINW_CH#{adc.to_s}"
    energy_width=event_hdu.header(energy_width_header).to_f
    unix_time_start=event_hdu["unixTime"][0].to_f
    unix_time_last=event_hdu["unixTime"][event_num].to_f
    observation_time=unix_time_last-unix_time_start

    for i in 0..event_num
      if adc_index[i].to_i==adc then
        energy=(energy_raw[i].to_f+(rand(-500..500).to_f/1000.0)*energy_width)/1000.0
        hist.Fill(energy)
      end
    end

    for i in 0..num-1
      scale=1.0/(observation_time*(bins[i+1]-bins[i]))
      bin_scaled=(hist.GetBinContent(i+1))*scale
      bin_scaled_error=(hist.GetBinError(i+1))*scale
      hist.SetBinContent(i+1, bin_scaled)
      hist.SetBinError(i+1, bin_scaled_error)
    end

    hist.SetTitle("")
    hist.GetXaxis().SetTitle("Eenergy (MeV)")
    hist.GetYaxis().SetTitle("Spectrum (counts s^{-1} MeV^{-1})")
    hist.SetStats(0)

    return hist
  end

  def ql_lightcurve(adc=0, min=0.3, max=30.0, width=1.0, name="ql_lightcurve")
    if @fits.length!=1 then
      puts "Error: Quick look mode can be used with a single fits file."
      exit 1
    end

    fits=@fits_ins[0]
    event_hdu=fits[0]["EVENTS"]
    event_num=event_hdu.getNRows()-1
    adc_index=event_hdu["boardIndexAndChannel"]
    energy_raw=event_hdu["energy"]
    energy_width_header="BINW_CH#{adc.to_s}"
    energy_width=event_hdu.header(energy_width_header).to_f
    unix_time=event_hdu["unixTime"]
    
    unix_time_start=event_hdu["unixTime"][0].to_f
    unix_time_last=event_hdu["unixTime"][event_num].to_f
    observation_time=unix_time_last-unix_time_start
    bin_num=(observation_time/width).ceil
    bin_start=unix_time_start
    bin_last=bin_start+width*bin_num.to_f

    hist=Root::TH1D.create(name, name, bin_num, bin_start, bin_last)

    for i in 0..eventNum
      if (adc==adc_index[i].to_i)&&(unix_time[i]>=bin_start)&&(unix_time[i]<=bin_last) then
        energy=(energy_raw[i]+(rand(-50..50).to_f/100.0)*energy_width)/1000.0
        if (energy>=min)&&(energy<max) then
          hist.Fill(unix_time[i].to_f)
        end
      end
    end

    scale_factor=1.0/width
    hist.Scale(scale_factor)
    hist.SetTitle("")
    hist.GetXaxis().SetTitle("Date (JST)")
    hist.GetYaxis().SetTitle("Count Rate (counts s^{-1})")
    hist.SetStats(0)
    hist.GetXaxis().SetTimeDisplay(1)
    hist.GetXaxis().SetTimeFormat("%H:%M")

    return hist
  end
  
  def lightcurve(adc=0, center, before, after, min=0.3, max=30.0, width=1.0, name="lightcurve")
    bin_before=(before/width).ceil
    bin_after=(after/width).ceil
    num=bin_before+bin_after
    range=[-1.0*width*bin_before.to_f, width*bin_after.to_f]

    all=Root::TH1F.create("all", "all", num, range[0], range[1])
    cnt=Root::TH1F.create("cnt", "cnt", num, range[0], range[1])
    raw=Root::TH1F.create("raw", "raw", num, range[0], range[1])
    rat=Root::TH1F.create("rat", "rat", num, range[0], range[1])
    cor=Root::TH1F.create(name, name, num, range[0], range[1])

    @fits.each do |fits_ins|

      event_hdu=fits_ins["EVENTS"]
      adc_index=event_hdu["boardIndexAndChannel"]
      energy_raw=event_hdu["energy"]
      unix_time=event_hdu["unixTime"]
      trigger_count=event_hdu["triggerCount"]
      energy_width_header="BINW_CH#{adc.to_s}"
      energy_width=event_hdu.header(energy_width_header).to_f
      event_num=event_hdu.getNRows()-1
      trigger_count_past=-1
      unix_time_past=0.0

      for i in 0..eventNum
        if (adc==adc_index[i].to_i)&&(unix_time[i].to_f-center>=range[0])&&(unix_time[i].to_f-center<=range[1]) then
          all.Fill(unix_time[i]-center)
          cnt.Fill(unix_time[i]-center)
          if triggerCountPast>=0 then
            delta_trigger_count=trigger_count[i].to_i-trigger_count_past-1
            if delta_trigger_count<0 then
              delta_trigger_count+=2**16
            end
            delta_trigger_count.times do
              probability=rand(0..10000).to_f/10000.0
              lost_time=unix_time_past+(unix_time[i].to_f-unix_time_past)*probability
              all.Fill(lost_time-center)
            end
          end
          trigger_count_past=trigger_count[i].to_i
          unix_time_past=unix_time[i].to_f
          energy=(energy_raw[i].to_f+(rand(-500..500).to_f/1000.0)*energy_width)/1000.0
          if (energy>=min)&&(energy<=max) then
            raw.Fill(unix_time[i].to_f-center)
          end
        end
      end
    end
    rat.Divide(all, cnt, 1.0, 1.0)
    cor.Multiply(raw, rat, 1.0, 1.0)
    cor.Scale(1.0/width)

    cor.SetTitle("")
    cor.GetXaxis().SetTitle("Time (second)")
    cor.GetYaxis().SetTitle("Count Rate (counts s^{-1})")
    cor.SertStats(0)
    return cor
  end

  def lc_unixtime(adc=0, center, before, after, min=0.3, max=30.0, width=1.0, name="lc_unixtime")
    bin_before=(before/width).ceil
    bin_after=(after/width).ceil
    num=bin_before+bin_after
    range=[center-1.0*width*bin_before.to_f, center+width*bin_after.to_f]

    all=Root::TH1F.create("all", "all", num, range[0], range[1])
    cnt=Root::TH1F.create("cnt", "cnt", num, range[0], range[1])
    raw=Root::TH1F.create("raw", "raw", num, range[0], range[1])
    rat=Root::TH1F.create("rat", "rat", num, range[0], range[1])
    cor=Root::TH1F.create(name, name, num, range[0], range[1])

    @fits.each do |fits_ins|

      event_hdu=fits_ins["EVENTS"]
      adc_index=event_hdu["boardIndexAndChannel"]
      energy_raw=event_hdu["energy"]
      unix_time=event_hdu["unixTime"]
      trigger_count=event_hdu["triggerCount"]
      energy_width_header="BINW_CH#{adc.to_s}"
      energy_width=event_hdu.header(energy_width_header).to_f
      event_num=event_hdu.getNRows()-1
      trigger_count_past=-1
      unix_time_past=0.0

      for i in 0..eventNum
        if (adc==adc_index[i].to_i)&&(unix_time[i].to_f>=range[0])&&(unix_time[i].to_f<=range[1]) then
          all.Fill(unix_time[i])
          cnt.Fill(unix_time[i])
          if triggerCountPast>=0 then
            delta_trigger_count=trigger_count[i].to_i-trigger_count_past-1
            if delta_trigger_count<0 then
              delta_trigger_count+=2**16
            end
            delta_trigger_count.times do
              probability=rand(0..10000).to_f/10000.0
              lost_time=unix_time_past+(unix_time[i].to_f-unix_time_past)*probability
              all.Fill(lost_time)
            end
          end
          trigger_count_past=trigger_count[i].to_i
          unix_time_past=unix_time[i].to_f
          energy=(energy_raw[i].to_f+(rand(-500..500).to_f/1000.0)*energy_width)/1000.0
          if (energy>=min)&&(energy<=max) then
            raw.Fill(unix_time[i].to_f)
          end
        end
      end
    end
    rat.Divide(all, cnt, 1.0, 1.0)
    cor.Multiply(raw, rat, 1.0, 1.0)
    cor.Scale(1.0/width)

    cor.SetTitle("")
    cor.GetXaxis().SetTitle("Time (second)")
    cor.GetYaxis().SetTitle("Count Rate (counts s^{-1})")
    cor.GetXaxis().SetTimeDisplay(1)
    cor.GetXaxis().SetTimeFormat("%H:%M")
    cor.SertStats(0)
    return cor
  end

  
  def spectrum(adc=0, start, duration, min=0.05, max=50.0, num=200, log=true, name="spectrum")

    range=[start, start+duration]

    if log==true then
      bins=define_bins(min, max, num)
    else
      bins=Root::DoubleArray.new(num+1)
      width=(max-min)/num.to_f
      for i in 0..num
        bins[i]=min+width*i.to_f
      end
      bins[0]=min
      bins[num]=max
    end
    hist=Root::TH1D.create(name, name, num, bins)

    @fits.each do |fits_ins|
      event_hdu=fits_ins[0]["EVENTS"]
      event_num=event_hdu.getNRows()-1
      adc_index=event_hdu["boardIndexAndChannel"]
      energy_raw=event_hdu["energy"]
      unix_time=event_hdu["unixTime"]
      trigger_count=event_hdu["triggerCount"]
      energy_width_header="BINW_CH#{adc.to_s}"
      energy_width=event_hdu.header(energy_width_header).to_f
      event_num=event_hdu.getNRows()-1
      trigger_count_past=-1
      unix_time_past=0.0
      all=0
      cnt=0
      
      for i in 0..event_num
        if (adc_index[i].to_i==adc)&&(unix_time[i].to_f>=range[0])&&(unix_time[i].to_f<=range[1]) then
          energy=(energy_raw[i].to_f+(rand(-500..500).to_f/1000.0)*energy_width)/1000.0
          hist.Fill(energy)
          cnt+=1
          all+=1
        end
        if adc_index[i].to_i==adc then
          if triggerCountPast>=0 then
            delta_trigger_count=trigger_count[i].to_i-trigger_count_past-1
            if delta_trigger_count<0 then
              delta_trigger_count+=2**16
            end
            delta_trigger_count.times do
              probability=rand(0..10000).to_f/10000.0
              lost_time=unix_time_past+(unix_time[i].to_f-unix_time_past)*probability
              if (lost_time>=range[0])&&(lost_time<=range[1]) then
                all+=1
              end
            end
          end
          trigger_count_past=trigger_count[i].to_i
          unix_time_past=unix_time[i].to_f
        end
      end
    end

    ratio=all.to_f/cnt.to_f
    for i in 0..num-1
      scale=ratio/(duration*(bins[i+1]-bins[i]))
      bin_scaled=(hist.GetBinContent(i+1))*scale
      bin_scaled_error=(hist.GetBinError(i+1))*scale
      hist.SetBinContent(i+1, bin_scaled)
      hist.SetBinError(i+1, bin_scaled_error)
    end

    hist.SetTitle("")
    hist.GetXaxis().SetTitle("Eenergy (MeV)")
    hist.GetYaxis().SetTitle("Spectrum (counts s^{-1} MeV^{-1})")
    hist.SetStats(0)

    return hist
    
  end

  def spec_src(adc=0, start, duration, bgd_start, bgd_duration, min=0.05, max=50.0, num=200, log=true, name="spec_src")

    range_src=[start, start+duration]
    range_bgd=[bgd_start, bgd_start+bgd_duration]
    
    if log==true then
      bins=define_bins(min, max, num)
    else
      bins=Root::DoubleArray.new(num+1)
      width=(max-min)/num.to_f
      for i in 0..num
        bins[i]=min+width*i.to_f
      end
      bins[0]=min
      bins[num]=max
    end
    src=Root::TH1D.create("src", "src", num, bins)
    bgd=Root::TH1D.create("bgd", "bgd", num, bins)
    hist=Root::TH1D.create(name, name, num, bins)

    @fits.each do |fits_ins|
      event_hdu=fits_ins[0]["EVENTS"]
      event_num=event_hdu.getNRows()-1
      adc_index=event_hdu["boardIndexAndChannel"]
      energy_raw=event_hdu["energy"]
      unix_time=event_hdu["unixTime"]
      trigger_count=event_hdu["triggerCount"]
      energy_width_header="BINW_CH#{adc.to_s}"
      energy_width=event_hdu.header(energy_width_header).to_f
      event_num=event_hdu.getNRows()-1
      trigger_count_past=-1
      unix_time_past=0.0
      all_src=0
      cnt_src=0
      all_bgd=0
      cnt_bgd=0
      
      for i in 0..event_num
        if (adc_index[i].to_i==adc)&&(unix_time[i].to_f>=range_src[0])&&(unix_time[i].to_f<=range_src[1]) then
          energy=(energy_raw[i].to_f+(rand(-500..500).to_f/1000.0)*energy_width)/1000.0
          src.Fill(energy)
          cnt_src+=1
          all_src+=1
        end
        if (adc_index[i].to_i==adc)&&(unix_time[i].to_f>=range_bgd[0])&&(unix_time[i].to_f<=range_bgd[1]) then
          energy=(energy_raw[i].to_f+(rand(-500..500).to_f/1000.0)*energy_width)/1000.0
          bgd.Fill(energy)
          cnt_bgd+=1
          all_bgd+=1
        end

        if adc_index[i].to_i==adc then
          if triggerCountPast>=0 then
            delta_trigger_count=trigger_count[i].to_i-trigger_count_past-1
            if delta_trigger_count<0 then
              delta_trigger_count+=2**16
            end
            delta_trigger_count.times do
              probability=rand(0..10000).to_f/10000.0
              lost_time=unix_time_past+(unix_time[i].to_f-unix_time_past)*probability
              if (lost_time>=range_src[0])&&(lost_time<=range_src[1]) then
                all_src+=1
              end
              if (lost_time>=range_bgd[0])&&(lost_time<=range_bgd[1]) then
                all_bgd+=1
              end
            end
          end
          trigger_count_past=trigger_count[i].to_i
          unix_time_past=unix_time[i].to_f
        end
      end
    end

    ratio_src=all_src.to_f/cnt_src.to_f
    ratio_bgd=all_bgd.to_f/cnt_bgd.to_f
    for i in 0..num-1
      scale_src=ratio_src/(duration*(bins[i+1]-bins[i]))
      bin_scaled_src=(src.GetBinContent(i+1))*scale_src
      bin_scaled_error_src=(src.GetBinError(i+1))*scale_src
      src.SetBinContent(i+1, bin_scaled_src)
      src.SetBinError(i+1, bin_scaled_error_src)
      scale_bgd=ratio_bgd/(bgd_duration*(bins[i+1]-bins[i]))
      bin_scaled_bgd=(sbgd.GetBinContent(i+1))*scale_bgd
      bin_scaled_error_bgd=(bgd.GetBinError(i+1))*scale_bgd
      bgd.SetBinContent(i+1, bin_scaled_bgd)
      bgd.SetBinError(i+1, bin_scaled_error_bgd)
    end
    hist.Add(src, bgd, 1.0, -1.0)
    hist.SetTitle("")
    hist.GetXaxis().SetTitle("Eenergy (MeV)")
    hist.GetYaxis().SetTitle("Spectrum (counts s^{-1} MeV^{-1})")
    hist.SetStats(0)

    return hist
    
  end

  def scatter(adc=0, center, before, after)
    range=[-1.0*before, after]
    time=Array.new
    min=Array.new
    max=Array.new
    
    @fits.each do |fits_ins|
      event_hdu=fits_ins[0]["EVENTS"]
      event_num=event_hdu.getNRows()-1
      adc_index=event_hdu["boardIndexAndChannel"]
      unix_time=event_hdu["unixTime"]
      precise=event_hdu["preciseTime"]
      pha_max=event_hdu["phaMax"]
      pha_min=event_hdu["phaMin"]
      for i in 0..eventNum
        unix_time_now=unix_time[i].to_f.floor+precise[i].to_f
        if (adc==adc_index[i].to_i)&&(unix_time_now-center>=range[0])&&(unix_time_now-center<=range[1]) then
          time << unix_time_now-center
          max << 5.0*(pha_max[i].to_f-2048.0)/2048.0
          min << 5.0*(pha_min[i].to_f-2048.0)/2048.0
        end
      end
    end
    
    graph=Array.new
    graph[0]=Root::TGraph.create(time, max)
    graph[1]=Root::TGraph.create(time, min)
    graph[0].SetTitle("")
    graph[0].GetXaxis().SetTitle("Time (second)")
    graph[0].GetYaxis().SetTitle("Voltage (V)")
    graph[1].SetTitle("")
    graph[1].GetXaxis().SetTitle("Time (second)")
    graph[1].GetYaxis().SetTitle("Voltage (V)")
    return graph
  end
end

