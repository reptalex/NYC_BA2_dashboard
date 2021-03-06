NYC BA.2 Comparative Epi Forecasts
================
Alex Washburne
April 10, 2022

## Comparing the New York City and United Kingdom BA.2 Outbreaks

This Github respository contains code to compare the ongoing New York City BA.2 outbreak with the preceding BA.2 outbreak in the United Kingdom, using [the same methods that proved effective at comparing Delta and Omicron outbreaks across regions](https://www.medrxiv.org/content/10.1101/2022.01.14.22269288v1) The idea is simple. We estimate the exponential growth rates of cases, `r(t)`, and carefully determine the start date, `t0` of an outbreak caused by the same variant of concern. We plot case growth rates as a function of "outbreak time" `t-t0`. The Delta and Omicron outbreaks followed similar arcs across regions for each variant of concern (VOC): case growth rates rose to a characteristic height and then fell to their peaks after a characteristic outbreak duration for each VOC.

Our forecasts for the NYC BA.2 outbreak are designed by the hypothesis that at any point in time our best guess for the future trends in case growth rates as a function of outbreak time for NYC will be the already-observed trends in case growth rates from the UK. Below is a gif showing how the UK and NYC outbreaks have evolved since late February, 2022.

![](README_files/figure-markdown_github/Outbreak%20Trajectory-1.gif)

To forecast cases, we extend the NYC outbreak from its latest observed value along a slope parallel to the growth rates observed in the UK, as shown in the dotted line below.

![](README_files/figure-markdown_github/Extrapolation%20of%20NYC%20using%20UK%20outbreak-1.png)

With predicted growth rates derived from extending the NYC trajectory parallel to the UK trajectory, we integrate forward in time from our current expected daily incidence, and convert the daily incidence forecasts into a the 7-day moving average forecast (dashed line with a 95% confidence intevral ribbon). Our forecast will look ahead 11 days from the last observed data point in NYC and include the probability of the NYC outbreak exceeding or dropping below the nearest thresholds for [the CDC's COVID-19 community levels](https://www.cdc.gov/coronavirus/2019-ncov/science/community-levels.html)

![](README_files/figure-markdown_github/Forecast-1.png)

This forecasting dashboard is designed as a regularly update of outbreak trajectories, with a transparent link between our forecasts and the underlying hypothesis that VOC-driven outbreaks have common characteristics. Like any hypothesis, ours can be rejected, and we provide `r(t)` trajectories to allow visual comparisons of observed outbreak features in relation to our hypothesis.

For more info, please see[the document here](https://github.com/reptalex/NYC_BA2_dashboard/tree/main/docs). If you have any questions or comments, please post here and I'll try to respond when time permits.
