# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require "physics_calculator/version"

Gem::Specification.new do |s|
  s.name        = "physics_calculator"
  s.version     = PhysicsCalculator::VERSION
  s.date        = '2013-09-04'
  s.authors     = ["Cristian Derr"]
  s.email       = ["cristian.j.derr@gmail.com"]
  s.license     = 'MIT'
  s.summary     = "Physics Calculator"
  s.description = "Gem that contains many useful methods using physics formulas. Many of them use the vector standard library. Formulas range from classical mechanics, quantum mechanics, electricity and magnetism, and special relativity."
  s.files         = `git ls-files`.split("\n")
  s.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  s.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  s.require_paths = ["lib"]

  # specify any dependencies here; for example:
  s.add_development_dependency "rspec"
  # s.add_runtime_dependency "rest-client"
end
