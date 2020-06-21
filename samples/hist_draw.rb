#!/usr/bin/env ruby
require "GROWTH"
require "RubyROOT"
include Root
include RootApp

fitsFile=["20170113_045008.fits.gz", "20170113_052012.fits.gz"]
growth=Growth.new(fitsFile)
center=growth.unixtime("2017-01-13 05:05:20")
spec=growth.spec_src(0, center-20.0, 20.0, center-180.0, 120.0, 0.3, 15.0, 50, true, "spec")

c0=Root::TCanvas.create("c0", "c0", 640, 480)
spec.Draw("e")
c0.SetLogx()
c0.SetLogy()
c0.Update()
run_app()
