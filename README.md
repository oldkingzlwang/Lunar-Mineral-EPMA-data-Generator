# Lunar-Mineral-EPMA-data-Generator
A useful tool to replicate realistic EPMA data of lunar minerals using computer.

## Background

Typical Electron Probe Micro-Analyzer (EPMA) experiments report compositional results for minerals as average values for each analytical spot, accompanied by the percentage error (Error%) for each element. Accounting for the analytical error introduced by EPMA is critical in Earth and planetary sciences, particularly when using geo-thermobarometers to estimate the formation or equilibrium pressure and temperature (P-T) conditions of rocks. Even small compositional variations can introduce significant uncertainties in P-T estimates, particularly for rocks from low-gravity celestial bodies (e.g., the Moon), as demonstrated by Wieser et al. (2023).

Recent geo-thermobarometer programs now integrate EPMA analytical errors (e.g., Ágreda-López et al., 2024; Wieser et al., 2022). A common and straightforward method assumes that the analytical error of each element is independent and follows a Gaussian normal distribution. This assumption is typically paired with Monte Carlo simulations, which generate random EPMA data by using the average values of the analytical spots and their 1σ errors as reported by the EPMA instrument.

However, this approach has a significant limitation: not all elements in natural minerals are independent of each other. For instance, FeO and MgO often exhibit strong negative correlations in natural minerals, and the elemental composition of a mineral must conform to its chemical formula. For example, in olivine, the sum of Fe, Mg, Ca, and Mn molar values should be approximately twice the molar value of Si. Ignoring these inter-element correlations can lead to the generation of unrealistic data that fails to represent natural mineral compositions.

To address this issue, I developed two MATLAB programs: the **Lunar Highland Mineral EPMA-data Generator (LHMEG)** and the **Lunar Mare Mineral EPMA-data Generator (LMMEG)**. These programs are designed to generate a specified number of random EPMA datasets while accounting for the mineral type, average major element compositions (including SiO<sub>2</sub>, TiO<sub>2</sub>, Al<sub>2</sub>O<sub>3</sub>, Cr<sub>2</sub>O<sub>3</sub>, FeO, MnO, MgO, CaO, Na<sub>2</sub>O, K<sub>2</sub>O, and P<sub>2</sub>O<sub>5</sub>), and their standard deviations. LHMEG is tailored for lunar highland rocks (including ferroan/magnesian anorthosites, magnesian-suite, alkali-suite), while LMMEG is designed for mare rocks (including mare basalts and picritic deposits). These programs ensure that the generated EPMA data respect inter-element correlations and adhere to the chemical formula of minerals, producing more realistic and mineral-specific EPMA datasets.

## Implementation of the program

This program is designed to model EPMA data for seven minerals in lunar rocks: pyroxene, olivine, plagioclase, ilmenite, Cr-rich spinel (Al<sub>2</sub>O<sub>3</sub> < 20 wt%), Mg-rich spinel (Al<sub>2</sub>O<sub>3</sub> > 40 wt%), and K-feldspar. Pyroxene is further subdivided into low-Ca pyroxene (CaO < 10 wt%) and high-Ca pyroxene (CaO > 10 wt%). To specify the mineral type, the ```mineral``` variable should be defined in either the **LHMEG** or **LMMEG** program.

To generate random EPMA compositions for these minerals, the program first requires predefined standard chemical characteristics. A preset data set of real EPMA measurements, termed "train data" (stored in the variable ```train_data```), is used to describe the chemical characteristics of minerals. The program calculates the Spearman correlation coefficients (```rho```) between every two elemental dimensions of this ```train_data```. The next step is to define the ```test_data```, which represents the desired “compositional shape” for the mineral to be modeled. This compositional shape is expressed as the mean (```test_mean```) and standard deviation (```test_std```) values of elemental compositions, derived from the ```test_data``` variable or inputted in the program manually.

The program allows for three different types of ```test_data``` input, controlled by the ```mode``` variable:

* **Mode = 0**: ```test_data``` is generated from a real EPMA database (e.g., ```Highlands.xlsx``` or ```Mare.xlsx```). The ```cvpartition``` function randomly splits the data into ```train_data``` and ```test_data``` according to a given proportion.
* **Mode = 1**: ```test_data``` is imported from a pre-defined file (```Test.xlsx```) containing EPMA compositional data, which the program will attempt to replicate.
* **Mode = 2**: The user manually provides the ```test_mean``` and ```test_std``` values (mean and standard deviations) without an explicit dataset.

Once the ```train_data``` and ```test_data``` are defined, the program calculates the Spearman correlation coefficients (```rho```) of the ```train_data```. The Spearman correlation coefficient is preferred over the Pearson coefficient because it does not require the data to be normally distributed.

Next, the program constructs the ```cov_matrix``` (covariance matrix) by multiplying ```rho``` with ```test_std```. The ```mvnrnd``` function is then used to generate random samples based on ```test_mean``` and ```cov_matrix```. The user specifies the number of desired EPMA samples (```num_samples```), which are then generated.

These generated samples are filtered to ensure they meet EPMA standards. Two constraints are applied:

  1. All compositional dimensions must have non-negative values.
  2. The sum of all compositional dimensions must lie within the range [98.5%, 101.5%].

The filtered samples are stored in the ```filtered_samples``` variable. Following this, the program performs a **Mann-Whitney U (MWU) test** to determine if the generated EPMA data in ```filtered_samples``` are consistent with the real EPMA data (```test_data```). The ```MWU_test``` function in the program tests the following criteria:

  1. **Pyroxene**: The molar ratios Ti/Al, Mg/Fe, Al(VI) + Fe<sup>3+</sup> – Na, and Na+K+Ca in the generated data should show no significant difference from the real EPMA data.
  2. **Olivine**: The molar Mg/Fe ratio in the generated data should show no significant difference from the real EPMA data.
  3. **Plagioclase and K-feldspar**: The An# and molar total amount in the generated data should show no significant difference from the real EPMA data.
  4. **Ilmenite and Spinel**: The Fe<sup>3+</sup>/(Fe<sup>3+</sup>+Ti), Fe<sup>2+</sup>/(Fe<sup>2+</sup>+Mn+Mg), and Ti/(Ti+Al+Cr) ratios in the generated data should show no significant difference from the real EPMA data.

The ```MWU_test``` function calls a series of sub-functions named "```test+mineral name```" (e.g., testPyroxene), which calculate the chemical formulas of the minerals using methods adapted from the **MinPlot** program (Walters, 2021). The program will then print the P-value for each test. If P-value > 0.05, the statement will indicate that there is no significant difference between the two datasets, meaning the generated data in ```filtered_samples``` accurately reflect the chemical nature of the real EPMA data and can be used.

If the ```test_data``` is predefined (i.e., when ```mode = 0``` or ```mode = 1```), the program will automatically generate comparative plots of the compositional dimensions using the ```gplotmatrix``` function. These plots display the real EPMA data (blue spots) and the generated data (red spots), allowing the user to visually assess whether the two datasets exhibit systematic differences.

Once all steps are complete, the program generates a results table and saves the generated EPMA data along with their molar compositions in a .xlsx file named “```Output + mineral name.xlsx```” in the working directory.

## How to Use the Program

  1. **Install Required Toolbox**: Ensure that the Statistics and Machine Learning Toolbox is successfully installed in MATLAB.
  2. **Select the Program**: For modeling minerals from lunar highlands, open ```LHMEG.m```. For modeling minerals from lunar mare, open ```LMMEG.m```.
  3. **Set Program Parameters**: Modify the following variables in the script to suit your needs:
     * ```mode```: Controls the input mode for generating EPMA data (explained below).
     * ```mineral```: Defines the mineral type (e.g., pyroxene, olivine).
     * ```numValidSamples```: Specifies the number of valid EPMA samples to be generated.
  4. **Input Data**:
     * If ```mode = 1```: Open the “```Test.xlsx```” file and replace the data with your own dataset.
     * If ```mode = 2```: Manually specify the mean (```test_mean```) and standard deviation (```test_std```) values for each element in lines 119-120 of the script.
  5. **Run the Program**: After setting up the program, click "Run" to execute the script and generate the modeled EPMA data.
  6. **Customization**: You can further customize the program by editing the code or the ```MWU_test``` function to adjust the testing criteria. You may also modify the contents of the two training databases (```Highland.xlsx``` and ```Mare.xlsx```) to include additional data or remove unnecessary entries, thus enhancing the model's flexibility.

## Remarks

This program is currently in its early stages of development and requires further refinement. Additionally, the mineral database should be expanded to cover a broader range of compositions for more accurate simulations.

If you would like to contribute to the program's development, or if you have any questions or feedback, please contact the author via email at ```zilong.wang@pku.edu.cn```.

## References
  * Ágreda-López M., Parodi V., Musu A., et al. 2024. Enhancing machine learning thermobarometry for clinopyroxene-bearing magmas. *Computers & Geosciences* 193: 105707. doi: ```10.1016/j.cageo.2024.105707```.
  * Walters J. B. 2022. MinPlot: A mineral formula recalculation and plotting program for electron probe microanalysis. *Mineralogia* 53: 51-66. doi: ```10.2478/mipo-2022-0005```.
  * Wieser P. E., Petrelli M., Lubbers J., et al. 2022. Thermobar: an open-source Python3 tool for thermobarometry and hygrometry. *Volcanica* 5(2): 349-384. doi: ```10.30909/vol.05.02.349384```.
  * Wieser P. E., Till C., Kent A., Gleeson M. 2023. Comment on ‘The magmatic architecture and evolution of the Chang’e-5 lunar basalts’. Preprint submitted to *EarthArxiv*. doi: ```10.31223/X5MM3B```.

