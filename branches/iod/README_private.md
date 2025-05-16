# Generating datasets for IOD branch

Make sure that the folder **`ODJobs_Simulated_Data_20001_20005_v2`** and the file **`20001_tim_tle`** are located at the same level as this readme. 

run `scripts/simulate_iod.py` to:

- Take in all of the ODJobs file names
- Reduce to 445 files (that we have truth info for)
- For each file, save t, ra, dec, px, py, pz as .csv files under `public/SAT_[XXXX].csv`
- Save the 445 truth satellite ID's under `private/sat_ids.csv`
- Import the truth data and convert the TLE's into (r, v) and t for each satellite using SSA code
- Saves the truth data from the first 5 satellites to `public/truth5.fits`
- Saves the complete truth information to `private/truth.fits`
- Saves the sample submission file for the first 5 satellites to `public/sample_submission5.yaml'
- Saves the sample submission for all 445 satellites to `private/sample_submission.yaml`

# ***Scores***

The `score_history.txt` file has the following format by column:

1) Competitor name
2) Name of sample file
3) Date
4) RMSE x
5) RMSE y
6) RMSE z
7) RMSE vx
8) RMSE vy
9) RMSE vz
10) % Predictions within FOV After 5 mins
11) % Predictions within FOV After 15 mins
12) % Predictions within FOV After 30 mins
13) % Predictions within FOV After 60 mins
14) CVM for accuracy after 5 min propagation
15) CVM for accuracy after 15 min propagation
16) CVM for accuracy after 30 min propagation
17) CVM for accuracy after 60 min propagation
18) CVM for accuracy of state estimate
19) Z score [0%, 20%]
20) Z score [20%, 40%]
21) Z score [40%, 60%]
22) Z score [60%, 80%]
23) Z score [80%, 100%]