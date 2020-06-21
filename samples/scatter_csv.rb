#!/usr/bin/env ruby
require "GROWTH"
require "RubyROOT"
include Root

fitsFile=["20170113_045008.fits.gz", "20170113_052012.fits.gz"]
growth=Growth.new(fitsFile)

center=growth.unixtime("2017-01-13 05:05:20")
scat=growth.scatter(0, center, 1.0, 1.0)
  
outputFile="scatter.csv"
num=scat[2].length
File.open(outputFile, "w") do |output|
  output.puts("x,phaMax,phaMin")
  for i in 0..num-1
    output.puts("#{scat[2][i]},#{scat[3][i]},#{scat[4][i]}")
  end
end
