using System;
using System.Collections.Generic;

namespace OMGI.AbsoluteReturn.TimeSeries
{

    public static class QuadraticVariation
    {
        /// <summary>
        /// Calculates realized quadratic variation of a series of prices/rates.
        /// </summary>
        /// <param name="dates"></param> // In ascending order
        /// <param name="data"></param> // ditto
        /// <param name="isPrice"></param>
        /// <param name="qv"></param>
        /// <param name="qvc"></param>
        /// <param name="qvNoisy"></param>
        /// <param name="iGamma">Estimate of the the variance of the observation error term, the noise term.</param>
        /// <param name="v"></param>
        /// <param name="deltas"></param>
        /// <param name="zBar"></param>
        public static void QuadVar(
            double[] timeDeltas,
            double[] observationDeltas,
            TuningParameters tuningParameters,
            ref double qv,     
            ref double zBar,
            ref double zHat)
        {
            int k = observationDeltas.Length;
            CalculateLocallyAveragedPriceDeltas(observationDeltas, k, out zBar, out zHat);

            qv = (zBar * zBar) - (0.5 * zHat);
            qv /= (k * tuningParameters.Lambda);
        }
        

        /// <summary>
        /// Calculate realized quadratic covariation of two series 
        /// </summary>
        public static void QuadCoVar(
            double[] observationDeltas1,
            double[] observationDeltas2,
            TuningParameters tuningParameters,
            double qv1,
            double zBar1,
            double v1,
            double qv2,
            double zBar2,
            double v2,
            ref double zHat,
            ref double qcov,
            ref double qcovc)
        {
            int k = observationDeltas1.Length;
            CalculateLocallyAveragedPriceDeltas(
                        observationDeltas1,
                        observationDeltas2,
                        k,
                        out zHat);


            qcov = (zBar1 * zBar2) - (0.5 * zHat);
            qcov /= (k * tuningParameters.Lambda); 

            qcovc = (IsJump(zBar1, v1) || IsJump(zBar2, v2)) ? 0.0 : qcov;
            qcovc /= (k * tuningParameters.Lambda); 
        }


        public static void CalculateCorrelations(
            double qv1,
            double qv2,
            double qcov,
            double qcovc,
            ref double corr,
            ref double corr_c,
            ref double corr_d,
            TuningParameters tuningParameters)
        {
            double n = Math.Sqrt(qv1 * qv2);
            corr = qcov / n;
            corr_c = qcovc / n;
            corr_d = (qcov - qcovc) / n;
        }



        /// <summary>
        /// Return the difference between two dates, expressed as the decimal proportion of a day, and with single-second precision.
        /// e.g. 24 hours = 1.0, 1 hour = 1 / 24, 1 second = 1 / (60 * 60 * 24)
        /// </summary>
        /// <param name="startTime"></param>
        /// <param name="endTime"></param>
        /// <returns></returns>
        public static double CalculateTimeDelta(DateTimeOffset startTime, DateTimeOffset endTime)
        {
            return Math.Sqrt((endTime - startTime).TotalSeconds / TimeSpan.FromDays(1).TotalSeconds);
        }


        /// <summary>
        /// Calculates the returns in the observed value between 2 observations.
        /// </summary>
        /// <param name="currentObservation"></param>
        /// <param name="priorObservation"></param>
        /// <param name="isPrice">True if the observations are prices, false if they are rates.</param>
        /// <returns></returns>
        public static double CalculateObservationDelta(double currentObservation, double priorObservation, bool isPrice)
        {
            if (isPrice)
            {
                return Math.Log(currentObservation) - Math.Log(priorObservation);
            }
            else
            {
                return currentObservation - priorObservation;
            }
        }


        /// <summary>
        /// Calculates locally averaged increments and locally averaged squared increments
        /// </summary>
        /// <param name="delta"></param>
        /// <param name="k">The number of observations over which to perform the averaging.</param>
        /// <param name="zBar">Locally averaged increments.</param>
        /// <param name="zHat">Locally averaged squared increments. The bias correction for the noise, it estimates its variance.</param>
        private static void CalculateLocallyAveragedPriceDeltas(
            IReadOnlyList<double> delta,
            int k,
            out double zBar,
            out double zHat)
        {
            zBar = 0.0;
            zHat = 0.0;
            int i = delta.Count - k + 2 - 1;
            for (int j = 1; j < k; j++)
            {
                double weight = g(j / (double)k); // if this is called "weight", what should we call g((j + 1) / (double)k)
                zBar += weight * delta[i + j - 1];
                zHat += Math.Pow((g((j + 1) / (double)k) - weight), 2) * delta[i + j - 1] * delta[i + j - 1];
            }
        }


        /// <summary>
        /// Scope to merge these 2 similarly named methods as code is duplicated
        /// </summary>
        /// <param name="delta1"></param>
        /// <param name="delta2"></param>
        /// <param name="k"></param>
        /// <param name="zHat"></param>
        private static void CalculateLocallyAveragedPriceDeltas(
            IReadOnlyList<double> delta1,
            IReadOnlyList<double> delta2,
            int k,
            out double zHat)
        {
            zHat = 0.0;
            int i = delta1.Count - k + 2 - 1;
            for (int j = 1; j < k; j++)
            {
                zHat += Math.Pow((g((j + 1) / (double)k) - g(j / (double)k)), 2) * delta1[i + j - 1] * delta2[i + j - 1];
            }
        }



        public static double CalculateTheta(double[] timeDeltas, int i, int k)
        {
            double theta = 0.0;
            int start = i - k + 1;
            for (int j = 0; j < k; j++)
            {
                int index = start + j;
                theta += Math.Sqrt(timeDeltas[index]);
            }
            return theta;
        }
        

        /// <summary>
        /// Localising function. 
        /// </summary>
        /// <param name="x">The value to localise.</param>
        /// <returns>A value that satifies g(0) = 0, g(1) = 0 and integral from 0, 1 of g^2 > 0.</returns>
        private static double g(double x)
        {
            return Math.Min(x, 1 - x);
        }


        /// <summary>
        /// Truncation function. Indicates if a value is likly to be a jump.
        /// </summary>
        /// <param name="x">Locally averaged price change.</param>
        /// <param name="sigmaContinuous">Standard deviation of the continuous part.</param>
        /// <returns>True if the price change is probably a jump, i.e. greater than the standard deviation, otherwise false.</returns>
        private static bool IsJump(double x, double sigmaContinuous)
        {
            return Math.Abs(x) > sigmaContinuous;
        }




        /// <summary>
        /// Normalises an accumalated quadratic variation of a bond future's price to basis point per day volatility
        /// </summary>
        /// <param name="interval">The time interval over which the <paramref name="qv"/>qv</param> has accumalated.
        /// <param name="qv">The accumalated quadratic variation.</param>
        /// <param name="ctdDuration">The modified duration of the future's cheapest to deliver bond.</param>
        /// <returns></returns>
        public static double Normalise(TimeSpan interval, double qv, double ctdDuration) // Specific to bond futures, not quad var in general, so does not belong in here
        {
            if(qv < 0.0)
            {
                return 0.0;
            }
            // The proportion of a day that is in our accumalation interval
            double intervalProportion = interval.TotalSeconds / TimeSpan.FromDays(1).TotalSeconds;
            // The rate at which quad variation is accumulating, per unit time
            double accumalationRate = qv / intervalProportion;
            // Take square root and multiply by duration of the CTD bond for future in question
            double normalisedQv = Math.Sqrt(accumalationRate) * ctdDuration;
            // Mutltiply by 100 to get a basis point per day volatility
            return normalisedQv * 100.0;
        }



        //private static void CalculateDecomposition(
        //    double[] deltas,
        //    double[] timeDeltas,
        //    double zBar,
        //    double zHat,
        //    int k,
        //    TuningParameters tuningParameters,
        //    double qv,
        //    ref double qvc,
        //    ref double qvNoisy,
        //    ref double iGamma,
        //    ref double v)
        //{

        //    if (k > 0)
        //    {
        //        //double u = CalculateTimeDelta(timeDeltas, i - k[i] + 1, k[i]);  
        //        //// double s0 = Statistics.StandardDeviation(new ArraySegment<double>(zBar, 0, i + 1));
        //        //stats.Push(zBar[i]); // should this go before the k test so zeroes values get pushed too?
        //        //double s0 = stats.StandardDeviation;
        //        //double v0 = tuningParameters.a0 * s0 * (Math.Pow(u, tuningParameters.w)); //%initial truncation parameter
        //        //double[] rets_trunc = new double[i + 1];
        //        //for (int j = 0; j < rets_trunc.Length; j++)
        //        //{
        //        //    rets_trunc[j] = Math.Min(deltas[j], v0);
        //        //}
        //        //double s = Statistics.StandardDeviation(rets_trunc);


        //        //v[i] = tuningParameters.a0 * s * (Math.Pow(u, tuningParameters.w)); // standard deviation of the continuous part
        //        //double minNonZeroAbs = deltas.Take(i + 1).Where(d => d != 0.0).MinimumAbsolute();
        //        //v[i] = Math.Max(v[i], double.IsNaN(minNonZeroAbs) ? v[i] : minNonZeroAbs);

        //        //qvNoisy[i] = Math.Pow(deltas[i], 2);
        //        //qv[i] = Math.Pow(zBar[i], 2) - 0.5 * zHat[i];
        //        //if (!IsJump(zBar[i], v[i]))
        //        //{
        //        //    qvc[i] = (Math.Pow(zBar[i], 2) - 0.5 * zHat[i]);
        //        //    iGamma[i] = zHat[i];
        //        //}
        //        //qv[i] = (1 / (k[i] * tuningParameters.Lambda)) * qv[i];
        //        //qvc[i] = (1 / (k[i] * tuningParameters.Lambda)) * qvc[i];
        //        //double theta = CalculateTheta(timeDeltas, i, k[i]);
        //        //iGamma[i] = (Math.Pow(theta, 2) / (2 * k[i] * tuningParameters.Lambda)) * iGamma[i];
        //    }
        //}
    }
}
