#!/usr/bin/env ruby
require "GROWTH"
require "RubyROOT"
include Root

fitsFile=["20170113_045008.fits.gz", "20170113_052012.fits.gz"]
growth=Growth.new(fitsFile)

center=growth.unixtime("2017-01-13 05:05:20")
scat=growth.scatter(0, center, 1.0, 1.0)
  
outputFile="scatter.root"
TFile.open(outputFile, "RECREATE") do |root|
  scat[0].Write("phaMax")
  scat[1].Write("phaMin")
end
