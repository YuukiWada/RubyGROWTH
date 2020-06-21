#!/usr/bin/env ruby
require "GROWTH"
require "RubyROOT"
include Root

fitsFile=["20170113_045008.fits.gz", "20170113_052012.fits.gz"]
growth=Growth.new(fitsFile)
center=growth.unixtime("2017-01-13 05:05:20")
spec=growth.spec_src(0, center-20.0, 20.0, center-180.0, 120.0, 0.3, 15.0, 50, true, "spec")

outputFile="spectrum.qdp"
File.open(outputFile, "w") do |output|
  output.puts("READ SERR 1 2")
  output.puts("@test.pco")
  output.puts("!")
  binNum=spec.GetNbinsX()
  for i in 1..binNum
    x=spec.GetBinCenter(i)
    ex=spec.GetBinWidth(i)/2.0
    y=spec.GetBinContent(i)
    ey=spec.GetBinError(i)
    output.puts("#{x} #{ex} #{y} #{ey}")
  end
end
