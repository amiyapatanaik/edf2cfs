// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_FFTW_H
#define SP_FFTW_H
#include <fftw3.h>

namespace sp
{
	///
	/// @defgroup fftw FFTW
	/// \brief One dimensional FFT functions using FFTW3 library.
	///
	/// \note If a single FFT is to be used the Armadillo version is faster. 
	/// FFTW takes longer time at the first calculation but is faster in the following loops 
	/// @{


	///
	/// \brief FFTW class.
	///
	/// Implements FFT functions for Armadillo types. For more info see [fftw.org](http://fftw.org/)
	/// 
	class FFTW
	{
	private:
		fftw_plan pl_fft;     ///< Real FFTW plan 
		fftw_plan pl_ifft;    ///< Real IFFTW plan 
		fftw_plan pl_fft_cx;  ///< Complex FFTW plan 
		fftw_plan pl_ifft_cx; ///< Complex IFFTW plan 
		int N;                ///< FFT length
		int alg;              ///< One of FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE, FFTW_WISDOM_ONLY see [FFTW plans](http://fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags)			
		int export_alg;       ///< Alg used for exporting wisdom
	public:
		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Constructor.
		/// @param _N FFT length
		/// @param _alg FFTW algorithm selection, Default FFTW_ESTIMATE
		////////////////////////////////////////////////////////////////////////////////////////////
		FFTW(int _N, int _alg = FFTW_ESTIMATE)
		{
			N = _N;
			alg = _alg;
			export_alg = FFTW_PATIENT;
			pl_fft = NULL;
			pl_ifft = NULL;
			pl_fft_cx = NULL;
			pl_ifft_cx = NULL;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Destructor.
		////////////////////////////////////////////////////////////////////////////////////////////
		~FFTW()
		{
			if (pl_fft != NULL) fftw_destroy_plan(pl_fft);
			if (pl_ifft != NULL) fftw_destroy_plan(pl_ifft);
			if (pl_fft_cx != NULL) fftw_destroy_plan(pl_fft_cx);
			if (pl_ifft_cx != NULL) fftw_destroy_plan(pl_ifft_cx);
			fftw_cleanup();
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief FFT of complex input.
		/// @param x Complex input data
		/// @param[out] Pxx Vector to hold complex FFT of length N 
		////////////////////////////////////////////////////////////////////////////////////////////
		void fft_cx(arma::cx_vec &x, arma::cx_vec &Pxx)
		{
			fftw_complex*  in = reinterpret_cast<fftw_complex*>(x.memptr());
			fftw_complex* out = reinterpret_cast<fftw_complex*>(Pxx.memptr());
			if (pl_fft_cx == NULL)
			{
				pl_fft_cx = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, alg);
				if (pl_fft_cx == NULL)
				{
					err_handler("Unable to create complex data FFTW plan");
				}
			}
			fftw_execute_dft(pl_fft_cx, in, out);
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief FFT of complex input.
		/// @returns Complex FFT of length N 
		/// @param x Complex input data
		////////////////////////////////////////////////////////////////////////////////////////////
		arma::cx_vec fft_cx(arma::cx_vec &x)
		{
			arma::cx_vec Pxx(N);
			fft_cx(x, Pxx);
			return Pxx;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Inverse FFT.
		/// @param Pxx Complex FFT
		/// @param[out] x Vector to hold complex data of length N 
		////////////////////////////////////////////////////////////////////////////////////////////
		void ifft_cx(arma::cx_vec &Pxx, arma::cx_vec &x)
		{
			fftw_complex*  in = reinterpret_cast<fftw_complex*>(Pxx.memptr());
			fftw_complex* out = reinterpret_cast<fftw_complex*>(x.memptr());
			if (pl_ifft_cx == NULL)
			{
				pl_ifft_cx = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, alg);
				if (pl_ifft_cx == NULL)
				{
					err_handler("Unable to create complex data IFFTW plan");
				}
			}
			fftw_execute_dft(pl_ifft_cx, in, out);
			x /= N;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Inverse FFT.
		/// @returns Complex data vector of length N 
		/// @param Pxx Complex FFT
		////////////////////////////////////////////////////////////////////////////////////////////
		arma::cx_vec ifft_cx(arma::cx_vec &Pxx)
		{
			arma::cx_vec x(N);
			ifft_cx(Pxx, x);
			return x;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief FFT of real input.
		/// @param x Input data
		/// @param[out] Pxx Vector to hold complex FFT of length N 
		////////////////////////////////////////////////////////////////////////////////////////////
		void fft(arma::vec &x, arma::cx_vec &Pxx)
		{
			double*        in = x.memptr();
			fftw_complex* out = reinterpret_cast<fftw_complex*>(Pxx.memptr());
			if (pl_fft == NULL)
			{
				pl_fft = fftw_plan_dft_r2c_1d(N, in, out, alg);
				if (pl_fft == NULL)
				{
					err_handler("Unable to create real data FFTW plan");
				}
			}

			fftw_execute_dft_r2c(pl_fft, in, out);
			int offset = static_cast<int>(ceil(N / 2.0));
			int n_elem = N - offset;
			for (int i = 0; i < n_elem; ++i)
			{
				Pxx(offset + i) = std::conj(Pxx(n_elem - i));
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief FFT of real input.
		/// @returns Complex FFT of length N 
		/// @param x Real input data
		arma::cx_vec fft(arma::vec &x)
		{
			arma::cx_vec Pxx(N);
			fft(x, Pxx);
			return Pxx;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Inverse FFT.
		/// @param Pxx Complex FFT
		/// @param[out] x Vector to hold real data of length N 
		////////////////////////////////////////////////////////////////////////////////////////////
		void ifft(arma::cx_vec &Pxx, arma::vec &x)
		{
			fftw_complex* in = reinterpret_cast<fftw_complex*>(Pxx.memptr());
			double*      out = x.memptr();
			if (pl_ifft == NULL)
			{
				pl_ifft = fftw_plan_dft_c2r_1d(N, in, out, alg);
				if (pl_ifft == NULL)
				{
					err_handler("Unable to create real data IFFTW plan");
				}
			}
			fftw_execute_dft_c2r(pl_ifft, in, out);
			x /= N;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Inverse FFT.
		/// @returns Real data vector of length N 
		/// @param Pxx Complex FFT
		////////////////////////////////////////////////////////////////////////////////////////////
		arma::vec ifft(arma::cx_vec &Pxx)
		{
			arma::vec x(N);
			ifft(Pxx, x);
			return x;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Import wisdom from string.
		/// @param wisd Wisdom string
		////////////////////////////////////////////////////////////////////////////////////////////
		void import_wisdom_string(const std::string wisd)
		{
			int res = fftw_import_wisdom_from_string(wisd.c_str());
			if (res == 0)
				err_handler("Unable to import wisdom from string!");
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Import wisdom from file.
		/// @param fname File name
		////////////////////////////////////////////////////////////////////////////////////////////
		void import_wisdom_file(const std::string fname)
		{
			int res = fftw_import_wisdom_from_filename(fname.c_str());
			if (res == 0)
				err_handler("Unable to import wisdom from file!");
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Export real FFT wisdom to file.
		/// @param fname File name
		////////////////////////////////////////////////////////////////////////////////////////////
		void export_wisdom_fft(const std::string fname)
		{
			fftw_plan pl_w = NULL;
			double *x_r = fftw_alloc_real(N);
			fftw_complex *x_cx1 = fftw_alloc_complex(N);

			// Replan using wisdom
			pl_w = fftw_plan_dft_r2c_1d(N, x_r, x_cx1, export_alg);
			if (pl_w == NULL)
			{
				err_handler("Unable to create real data FFTW plan");
			}

			// Export 
			if (fftw_export_wisdom_to_filename(fname.c_str()) == 0)
			{
				err_handler("Could not export wisdom to file!");
			}

			fftw_destroy_plan(pl_w);
			fftw_free(x_r);
			fftw_free(x_cx1);
		}
		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Export real IFFT wisdom to file.
		/// @param fname File name
		////////////////////////////////////////////////////////////////////////////////////////////
		void export_wisdom_ifft(const std::string fname)
		{
			fftw_plan pl_w = NULL;
			double *x_r = fftw_alloc_real(N);
			fftw_complex *x_cx1 = fftw_alloc_complex(N);

			// Replan using wisdom
			pl_w = fftw_plan_dft_c2r_1d(N, x_cx1, x_r, export_alg);
			if (pl_w == NULL)
			{
				err_handler("Unable to create real data FFTW plan");
			}

			// Export 
			if (fftw_export_wisdom_to_filename(fname.c_str()) == 0)
			{
				err_handler("Could not export wisdom to file!");
			}

			fftw_destroy_plan(pl_w);
			fftw_free(x_r);
			fftw_free(x_cx1);
		}
		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Export complex FFT wisdom to file.
		/// @param fname File name
		////////////////////////////////////////////////////////////////////////////////////////////
		void export_wisdom_fft_cx(const std::string fname)
		{
			fftw_plan pl_w = NULL;
			fftw_complex *x_cx1 = fftw_alloc_complex(N);
			fftw_complex *x_cx2 = fftw_alloc_complex(N);

			// Replan using wisdom
			pl_w = fftw_plan_dft_1d(N, x_cx1, x_cx2, FFTW_FORWARD, export_alg);
			if (pl_w == NULL)
			{
				err_handler("Unable to create complex data FFTW plan");
			}

			// Export 
			if (fftw_export_wisdom_to_filename(fname.c_str()) == 0)
			{
				err_handler("Could not export wisdom to file!");
			}

			fftw_destroy_plan(pl_w);
			fftw_free(x_cx1);
			fftw_free(x_cx2);
		}
		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Export complex IFFT wisdom to file.
		/// @param fname File name
		////////////////////////////////////////////////////////////////////////////////////////////
		void export_wisdom_ifft_cx(const std::string fname)
		{
			fftw_plan pl_w = NULL;
			fftw_complex *x_cx1 = fftw_alloc_complex(N);
			fftw_complex *x_cx2 = fftw_alloc_complex(N);

			// Replan using wisdom
			pl_w = fftw_plan_dft_1d(N, x_cx1, x_cx2, FFTW_BACKWARD, export_alg);
			if (pl_w == NULL)
			{
				err_handler("Unable to create complex data IFFTW plan");
			}

			// Export 
			if (fftw_export_wisdom_to_filename(fname.c_str()) == 0)
			{
				err_handler("Could not export wisdom to file!");
			}

			fftw_destroy_plan(pl_w);
			fftw_free(x_cx1);
			fftw_free(x_cx2);
		}
	};
	/// @}

} // end namespace
#endif