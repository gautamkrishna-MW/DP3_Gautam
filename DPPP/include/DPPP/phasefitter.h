/**
 * @file phasefitter.h Implements TEC model phase filter @ref PhaseFitter.
 * @author André Offringa
 * @date 2016-04-06
 */

#ifndef PHASE_FITTER_H
#define PHASE_FITTER_H

#include <vector>
#include <cmath>
#include <cstring>

/**
 * Phase fitter that can force phase solutions over frequency onto a TEC model.
 * To use:
 * - Construct and set the frequencies (with @ref FrequencyData() ) and if already possible
 *   the weights (@ref WeightData()).
 * - Perform a calibration iteration.
 * - Set the phase values (@ref PhaseData()) and, if not yet done, the weights (@ref WeightData())
 *   to the current phase solutions / weights of a single antenna.
 * - Call @ref FitDataToTECModel().
 * - Read the new phase values.
 * - When fitting multiple polarizations, the above steps should be done twice for each
 *   diagonal value of the Jones matrix. Also repeat for all antennae.
 * - Continue the iteration with the new phase values.
 */
class PhaseFitter
{
 public:
  /**
   * Construct a phase fitter for the given number of channels.
   * Weights are initialized to unity. The fitting accuracy is
   * initialized to 1e-6.
   * @param channelCount number of channels.
   */
  PhaseFitter(size_t channelCount) :
    _phases(channelCount, 0.0),
    _frequencies(channelCount, 0.0),
    _weights(channelCount, 1.0),
    _fittingAccuracy(1e-6)
  { }

  ~PhaseFitter() = default;
  PhaseFitter(const PhaseFitter&) = delete;
  PhaseFitter& operator=(const PhaseFitter&) = delete;

  /**
   * Fits the given phase values to a TEC model. This function is
   * robust even when phase wrapping occurs.
   * After this call, the @ref PhaseData() satisfy the TEC model.
   * The TEC model has two parameters and are fitted as described in
   * @ref FitTECModelParameters().
	 * @param alpha Found value for the alpha parameter.
	 * @param beta Found value for the beta parameter.
	 *
	 * @returns Cost of the found solution.
   */
  double FitDataToTECModel(double& alpha, double& beta);

  /**
   * Fits the given phase values to a TEC model using prior estimates of the 
	 * model parameters. This method is similar to @ref FitDataToTECModel(),
	 * except that it will use the provided initial values of alpha and
	 * beta to speed up the solution. When the provided initial values
	 * are not accurate, the fit might not be accurate.
	 * 
	 * @todo No fast method has been implemented -- instead it will perform a full parameter search.
	 * @param alpha Estimate of alpha parameter on input, found value on output.
	 * @param beta Estimate of beta parameter on input, found value on output.
	 *
	 * @returns Cost of the found solution.
   */
  double FitDataToTECModelWithInitialValues(double& alpha, double& beta)
	{
		return FitDataToTECModel(alpha, beta);
	}

  /**
   * Fit the data and get the best fitting parameters. The model
   * used is f(nu) = alpha/nu + beta, with possible 2pi wrappings in the
   * data.
   * The phase data is not changed.
   * The alpha parameter is linearly related to the TEC. The beta parameter
   * is a constant phase offset, given in radians. 
   * The fitting algorithm uses a combination of brute force searching and
   * binary-like searching (ternary search).
	 * @param alpha Will be set to the fitted value for the alpha parameter
	 * (value on input is not used).
	 * @param beta Will be set to the fitted value for the beta parameter
	 * (value on input is not used).
   */
  void FitTECModelParameters(double& alpha, double& beta) const;
  
	/**
	 * Get a pointer to the array of phase values. This array should be filled with the
	 * phases of solutions before calling one of the fit methods. @ref FitDataToTECModel()
	 * sets this array to the fitted phases.
	 * Normally, these values should be set to std::arg(z) or atan2(z.imag(), z.real()), where
	 * z is a complex solution for one polarizations. All phases should correspond to the same
	 * polarizations, i.e., different polarizations (xx/yy/ll/rr, etc.) should be independently
	 * fitted.
	 * @returns Array of @ref Size() doubles with the phases.
	 */
  double* PhaseData() { return _phases.data(); }
  
  /**
	 * Get a constant pointer to the array of values.
	 * @returns Constant array of @ref Size() doubles with the phases.
	 */
  const double* PhaseData() const { return _phases.data(); }

  /**
	 * Get a pointer to the array of frequency values. This array should be set to the
	 * frequency values of the channels, such that FrequencyData()[ch] is the frequency
	 * corresponding to the phase value PhaseData()[ch]. The fitter will not change this
	 * array.
	 * @returns Array of @ref Size() doubles with the frequencies in Hz.
	 */
  double* FrequencyData() { return _frequencies.data(); }
  
  /**
	 * Constant frequency data.
	 * @returns Constant array of @ref Size() doubles with the frequencies in Hz.
	 */
  const double* FrequencyData() const { return _frequencies.data(); }

  /**
   * This array should be filled with the weights
	 * of the channel phase solutions. If the solver supports weights during solving, this value should
	 * be set to the sum of weights of all visibilities that are used for the solution of
	 * this channel. If the solver does not support weights, it should be set to the number
	 * of unflagged visibilities used by the solver to generate the corresponding phase. Another
	 * way of saying this, is that the weights should be set to the esimated inverse variance of the phase
	 * solutions.
	 * 
	 * The use of weights will make sure that noisy channels do not bias the result. Weighting is
	 * for example helpful to avoid that the few remaining samples in a badly RFI contaminated
	 * channels cause the fit to be inaccurate.
	 * 
	 * While the weights could be different for each antenna solution, generally the weight of a channel
	 * is constant over the antennas. The latter implies that the weights can be set once before
	 * the solution starts, and only the @ref PhaseData() need to be changed within solution
	 * iterations.
	 * 
	 * The weights are initially set to one.
	 *
	 * @returns Array of @ref Size() doubles with the weights.
   */
  double* WeightData() { return _weights.data(); }
  
  /**
	 * Constant array of weights, as described above.
	 * @returns Constant array of @ref Size() doubles with weights.
	 */
  const double* WeightData() const { return _weights.data(); }

  /**
	 * Number of channels used for the fitter.
	 * @returns Number of channels.
	 */
  size_t Size() const { return _phases.size(); }

  /**
	 * Get the fitting accuracy. The fitter will stop once this accuracy is approximately reached.
	 * The default value is 1e-6.
	 * @returns Fitting accuracy.
	 */
  double FittingAccuracy() const { return _fittingAccuracy; }
  
  /**
	 * Change the fitting accuracy. See @ref FittingAccuracy().
	 * @param newAccuracy New accuracy.
	 */
  void SetFittingAccuracy(double newAccuracy) { _fittingAccuracy = newAccuracy; }

  /**
	 * Evaluate the cost function for given TEC model parameters. The higher the
	 * cost, the worser the data fit the given parameters.
	 * @param alpha TEC parameter
	 * @param beta Phase offset parameter
	 * @returns sum of | nu_i / alpha + beta - theta_i |, and each sum term is
	 * phase unwrapped.
	 */
	double TECModelCost(double alpha, double beta) const;
  
	/**
	 * Like @ref TECModelFunc(), but 2-pi wrapped.
	 * @param nu Frequency in Hz
	 * @param alpha TEC parameter (in undefined units)
	 * @param beta Phase offset parameter (in radians)
	 * @returns | nu_i / alpha + beta | % 2pi
	 */
	static double TECModelFunc(double nu, double alpha, double beta)
	{
		return alpha / nu + beta;
	}
  
	/**
	 * Like @ref TECModelFunc(), but 2-pi wrapped.
	 * @param nu Frequency in Hz
	 * @param alpha TEC parameter (in undefined units)
	 * @param beta Phase offset parameter (in radians)
	 * @returns | nu_i / alpha + beta | % 2pi
	 */
  static double TECModelFuncWrapped(double nu, double alpha, double beta)
  {
    return fmod(alpha / nu + beta, 2*M_PI);
  }
 private:
  std::vector<double> _phases, _frequencies, _weights;
  double _fittingAccuracy;
  
  double fitTECModelBeta(double alpha, double betaEstimate) const;
  void bruteForceSearchTECModel(double& lowerAlpha, double& upperAlpha, double& beta) const;
  double ternarySearchTECModelAlpha(double startAlpha, double endAlpha, double& beta) const;
  void fillDataWithTECModel(double alpha, double beta);
    
};

#endif
