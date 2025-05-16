## ***v2.0 updates:***
- Added "Summary" tab
- Added "Calibrate Sidereal" tab
- Changed various plot formats (added donut charts, added colroing by score to certain xy plots, etc)
- Added a new plot to IOD tab for the CVM (Cramer-von Mises statistic), to show the state estimate accuracy.
- Updated sample scores to better show off dashboard capabilities (while we wait on more real scores)

## ***v2.1 updates:***
- Added new plots to IOD tab: 
    - We now have 2 separate scores (with corresponding plots) for predicted accuracy downrange: scores using our propagator (which we call "standardized"), as well as scores using the competitor's propagator (if they have supplied us the optional requirements).

# ***Instructions:***

From this directory, run:

`bokeh serve --show dashboard`

optional arguments:

- `--dev` 

    - Including the --dev flag in the command makes Bokeh automatically reload 
    when any file changes are detected within the dashboard directory

# ***Summary Page Info:***

The summary page of the dashboard currently chooses one score from each branch to highlight the top scores from. Here are the current scores chosen by branch:

- IOD:
    - The percentage of satellites that are within the telescope FOV 30 minutes downrange.
- Detect Target:
    - The percentage of satellites correctly detected (by comparing determined location to truth data).
- Detect Sidereal:
    - The percentage of satellites correctly detected (by comparing determined location to truth data).
- Calibrate Sidereal:
    - The percentage of satellites correctly detected (by comparing determined location to truth data).