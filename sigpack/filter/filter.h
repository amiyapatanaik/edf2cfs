// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_FILTER_H
#define SP_FILTER_H
namespace sp
{
    ///
    /// @defgroup filter Filter
    /// \brief FIR/MA and IIR/ARMA filter functions.
    /// @{
        
    ///
    /// \brief FIR/MA filter class.
    ///
    /// Implements FIR/MA filter functions as \f[  y(n) = \sum_{k=0}^{N-1}{b_kx(n-k)}=b_0x(n)+b_1x(n-1)+...+b_{N-1}x(n-(N-1))\f]
    /// where N is the number of taps in the FIR filter. The filter order is N-1.  
    /// 
    template <class T1, class T2, class T3>
    class FIR_filt
    {
    private:
        int N;                ///< Nr of filter taps
        int cur_p;            ///< Pointer to current sample in buffer
        arma::Col<T1> buf;    ///< Signal buffer
        arma::Col<T2> b;      ///< Filter coefficients
    public:
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Constructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        FIR_filt(){}
        
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Destructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        ~FIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Clears the internal states and pointer.
        ////////////////////////////////////////////////////////////////////////////////////////////
        void clear(void)
        {
            buf.zeros();
            cur_p = 0;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Sets coefficients in FIR filter.
        /// The internal state and pointers are cleared
        /// @param _b Filter coefficients \f$ [b_0 ..b_{N-1}] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_coeffs(arma::Col<T2> &_b)
        {
            N = _b.size();
            buf.set_size(N);
            this->clear();
            b = _b;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Updates coefficients in FIR filter without clearing the internal states.
        /// @param _b Filter coefficients \f$ [b_0 ..b_{N-1}] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void update_coeffs(arma::Col<T2> &_b)
        {
            b = _b;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter operator.
        /// @returns Filtered output
        /// @param in Input sample
        ////////////////////////////////////////////////////////////////////////////////////////////
        T3 operator()(T1 & in)
        {
            T3 out=0;
            int p = 0;
            buf[cur_p] = in;                    // Insert new sample
            for(  int n = cur_p; n < N; n++)
                out += b[p++]*buf[n];           // Calc upper part
            for(int n = 0; n < cur_p; n++)
                out += b[p++]*buf[n];           // ... and lower

            cur_p--;                            // Move insertion point
            if (cur_p < 0) cur_p = N-1;

            return out;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter function.
        /// @returns Filtered output
        /// @param in Input vector
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Col<T3> filter(arma::Col<T1> & in)
        {
            long int sz = in.size();
            arma::Col<T3> out(sz);
            for(long int n=0;n<sz;n++)
                out[n] = this->operator()(in[n]);
            return out;
        }
    };


    ///
    /// \brief IIR/ARMA filter class.
    ///
    /// Implements IIR/ARMA filter functions as \f[  a_0y(n) = b_0x(n)+b_1x(n-1)+...+b_{N-1}x(n-(N-1))-a_1y(n-1)-...-a_{M-1}y(n-(M-1))\f]
    /// where N is the number of taps in the FIR filter part and M is the number of taps in the IIR filter. The filter order is (N-1,M-1)  
    /// 
    template <class T1, class T2, class T3>
    class IIR_filt
    {
    private:
        int N;                ///< Nr of MA filter taps
        int M;                ///< Nr of AR filter taps
        int b_cur_p;          ///< Pointer to current sample in MA buffer
        int a_cur_p;          ///< Pointer to current sample in AR buffer
        arma::Col<T2> b;      ///< MA Filter coefficients
        arma::Col<T2> a;      ///< AR Filter coefficients
        arma::Col<T1> b_buf;  ///< MA Signal buffer
        arma::Col<T1> a_buf;  ///< AR Signal buffer
    public:
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Constructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        IIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Destructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        ~IIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Clears the internal states and pointers.
        ////////////////////////////////////////////////////////////////////////////////////////////
        void clear(void)
        {
            b_buf.zeros();
            a_buf.zeros();
            b_cur_p = 0;
            a_cur_p = 0;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Sets coefficients in IIR filter.
        /// The internal state and pointers are cleared
        /// @param _b Filter coefficients \f$ [b_0 ..b_N] \f$
        /// @param _a Filter coefficients \f$ [a_0 ..a_M] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_coeffs(arma::Col<T2> &_b,arma::Col<T2> &_a)
        {
            N = _b.size();
            M = _a.size();
            b_buf.set_size(N);
            a_buf.set_size(M);
            this->clear();
            b = _b/_a[0];
            a = _a/_a[0];
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Updates coefficients in filter without clearing the internal states.
        /// @param _b Filter coefficients \f$ [b_0 ..b_N] \f$
        /// @param _a Filter coefficients \f$ [a_0 ..a_M] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void update_coeffs(arma::Col<T2> &_b,arma::Col<T2> &_a)
        {
            b = _b/_a[0];
            a = _a/_a[0];
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter operator.
        /// @returns Filtered output
        /// @param in Input sample
        ////////////////////////////////////////////////////////////////////////////////////////////
        T3 operator()(T1 & in)
        {
            T3 out=0;
            int p = 0;

            // MA part
            b_buf[b_cur_p] = in;                // Insert new sample
            for(int n = b_cur_p; n < N; n++)
                out += b[p++]*b_buf[n];         // Calc upper part
            for(int n = 0; n < b_cur_p; n++)
                out += b[p++]*b_buf[n];         // ... and lower

            b_cur_p--;                          // Move insertion point
            if (b_cur_p < 0) b_cur_p = N-1;

            // AR part
            p=1;
            for(int m = a_cur_p+1; m < M; m++)
                out -= a[p++]*a_buf[m];         // Calc upper part
            for(int m = 0; m < a_cur_p; m++)
                out -= a[p++]*a_buf[m];         // ... and lower

            a_buf[a_cur_p] = out;		        // Insert output

            a_cur_p--;
            if (a_cur_p < 0) a_cur_p = M-1;     // Move insertion point

            return out;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter function.
        /// @returns Filtered output
        /// @param in Input vector
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Col<T3> filter(arma::Col<T1> & in)
        {
            long int sz = in.size();
            arma::Col<T3> out(sz);
            for(long int n=0;n<sz;n++)
                out[n] = this->operator()(in[n]);
            return out;
        }
    };

    ///
    /// Filter design functions
    ///

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief FIR design functions.
    /// FIR design using windows method (hamming window).
    /// NB! Returns size N+1
    /// @returns b Filter coefficients \f$ [b_0 ..b_N] \f$
    /// @param N Filter order
    /// @param f0 Filter cutoff frequency in interval [0..1]
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec fir1(int N, double f0)
    {
        arma::vec b(N+1), h(N+1);
        h = hamming(N+1);
        double b_sum=0;
        for (int i=0;i<N+1;i++)
        {
            b[i] = h[i]*sinc(f0*(i-N/2.0));
            b_sum += b[i];
        }
        b = b/b_sum;
        return b;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Fractional delay function.
    /// Fractional delay filter design using windowed sinc method.
    /// Actual delay is N/2+fd samples for even nr of taps and
    /// (N-1)/2+fd for odd nr of taps
    /// Best performance if -1 < fd < 1
    /// @param N Filter length
    /// @param fd Fractional delay
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec fd_filter( const int N, double fd )
    {
        arma::vec h(N);
        arma::vec w = blackmanharris(N);
        if( N % 2 == 1 ) fd = fd-0.5; // Offset for odd nr of taps
        for(int n=0;n<N;n++)
        {
            h(n) = w(n)*sinc(n-N/2.0-fd);
        }
        h = h/sum(h);  // Normalize gain

        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Frequency response function.
    /// Calculates the frequency response
    /// @param b FIR/MA filter coefficients
    /// @param a IIR/AR filter coefficients
    /// @param M Number of evaluation points, Default 512
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::cx_vec freq( const arma::vec b, const arma::vec a, const int M=512)
    {
        arma::cx_vec h(M);
        int Nb = b.size();
        int Na = a.size();
        std::complex<double> b_tmp,a_tmp,i(0,1);
        for(int m=0;m<M;m++)
        {
            b_tmp=std::complex<double>(b(0),0);
            for(int nb=1;nb<Nb;nb++) 
                b_tmp+= b(nb)*(cos(nb*PI*m/M)-i*sin(nb*PI*m/M));
            a_tmp=std::complex<double>(a(0),0);
            for(int na=1;na<Na;na++) 
                a_tmp+= a(na)*(cos(na*PI*m/M)-i*sin(na*PI*m/M));
            h(m) = b_tmp/a_tmp;
        }
        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Frequency magnitude response function.
    /// Calculates the frequency magnitude response
    /// @param b FIR/MA filter coefficients
    /// @param a IIR/AR filter coefficients
    /// @param M Number of evaluation points, Default 512
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec freqz( const arma::vec b, const arma::vec a, const int M=512)
    {
        arma::cx_vec f = freq(b,a,M);
        return abs(f);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Frequency phase response function.
    /// Calculates the frequency phase response
    /// @param b FIR/MA filter coefficients
    /// @param a IIR/AR filter coefficients
    /// @param M Number of evaluation points, Default 512
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec phasez( const arma::vec b, const arma::vec a, const int M=512)
    {
        arma::cx_vec f = freq(b,a,M);
        return angle(f);
    }
    /// @}

} // end namespace
#endif