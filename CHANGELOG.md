# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed
Changed how ods files (netcdf version 5) were interacted with from the old hdfsd utility to the
new netcdf one.
### Fixed
The default settings only read 0 6 12 18
They now read every 3 hours

Default settings were failing for ozone/radiance.
Added a check to determine the types of obs being processed based up on the current directory and changed the options.dsynhhs variables in obs2html accordingly.

### Removed

### Deprecated

