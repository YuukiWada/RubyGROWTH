#!/usr/bin/env ruby
require "GROWTH"
require "RubyROOT"
include Root
include RootApp

fitsFile=["20170113_045008.fits.gz", "20170113_052012.fits.gz"]
growth=Growth.new(fitsFile)

center=growth.unixtime("2017-01-13 05:05:20")
scat=growth.scatter(0, center, 1.0, 1.0)
  
c0=Root::TCanvas.create("c0", "c0", 640, 480)
scat[0].GetYaxis().SetRangeUser(-5.0, 5.0)
scat[0].Draw("ap")
scat[1].Draw("p")
c0.Update()
run_app()
