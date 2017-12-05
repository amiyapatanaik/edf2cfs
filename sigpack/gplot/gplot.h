// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_GPLOT_H
#define SP_GPLOT_H
namespace sp
{
    ///
    /// @defgroup gplot GPlot
    /// \brief Collection of Gnuplot functions
    /// @{
        
    ///
    /// \brief Gnuplot class.
    ///
    /// Implements a class for streaming data to Gnuplot using a pipe.
    /// Inspiration from https://code.google.com/p/gnuplot-cpp/
    ///
    /// Verified with Gnuplot 4.6.5 for Win64 and Linux.
    /// \note In Windows only one class is allowed. Using multiple figures are controlled by a figure number. In Linux we may use one instance per figure.
    /// 
    class gplot
    {
    private:
        FILE           *gnucmd;          ///< File handle to pipe
        std::string    linestyle;        ///< Linestyle string as in Gnuplot doc.

    public:
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Constructor.
        ///
        /// Opens a pipe to gnuplot program. Make sure it is found/accessable by the system.
        ////////////////////////////////////////////////////////////////////////////////////////////
        gplot::gplot()
        {
        #if defined(WIN32)
            gnucmd = _popen("gnuplot -persist 2> NUL","w");
            #define GP_TERM "win"
        #elif defined(unix)
            //gnucmd = popen("gnuplot -persist &> /dev/null","w");
            gnucmd = popen("gnuplot -persist","w");
            #define GP_TERM "x11"
            //#elif defined(_APPLE_)
            //            gnucmd = popen("gnuplot -persist &> /dev/null","w");
            //#define GP_TERM "aqua"
        #else
            #error Only Windows and Linux/Unix is supported so far!
        #endif
            if(!gnucmd)
            {
                err_handler("Could not start gnuplot");
            }
            setvbuf(gnucmd, NULL, _IOLBF, 512);

            // Set global params
            linestyle = "lines";
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Destructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        gplot::~gplot()
        {
        #if defined(WIN32)
            _pclose(gnucmd);
        #elif defined(unix)
            pclose(gnucmd);
        #endif
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Send command to Gnuplot pipe.
        /// @param cmdstr  Command string
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::send2gp(const std::string &cmdstr)
        {
            std::string tmp=cmdstr+"\n";
            std::fputs(tmp.c_str(), gnucmd );
            //std::cout << tmp.c_str() << std::endl;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Sets the active figure.
        /// @param fig  Figure number
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::figure(const int fig)
        {
            std::ostringstream tmp_s;
            tmp_s << "set term " << GP_TERM << " " << fig;
            send2gp(tmp_s.str());
            send2gp("reset");
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Configure the figure used Windows environment.
        /// @param fig     Figure number
        /// @param name    Window name
        /// @param x       x position of upper left corner
        /// @param y       y position of upper left corner
        /// @param width   width of window
        /// @param height  height of window        
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::window(const int fig, const std::string &name,const int x,const int y,const int width,const int height)
        {
            std::ostringstream tmp_s;
            tmp_s << "set term " << GP_TERM << " " << fig << " title \"" << name << "\" position " << x << "," << y << " size " << width << "," << height;
            send2gp(tmp_s.str());
            send2gp("reset");
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Configure the figure/window - used in Linux environment where no figure numbers are needed.
        /// @param name    Window name
        /// @param x       x position of upper left corner
        /// @param y       y position of upper left corner
        /// @param width   width of window
        /// @param height  height of window        
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::window(const std::string &name,const int x,const int y,const int width,const int height)
        {
            window(0,name,x,y,width,height);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Set linestyle
        /// @param style  lines | points | linespoints | dots | steps         
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::set_linestyle(const std::string& style)
        {
            linestyle = style;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Set label for X-axis.
        /// @param label label string         
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::xlabel(const std::string& label)
        {
            std::ostringstream tmp_s;
            tmp_s << "set xlabel \"" << label << "\" ";
            send2gp(tmp_s.str());
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Set label for X-axis.
        /// @param label label string         
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::ylabel(const std::string& label)
        {
            std::ostringstream tmp_s;
            tmp_s << "set ylabel \"" << label << "\" ";
            send2gp(tmp_s.str());
        }

		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Set label at position x,y.
		/// @param x x value         
		/// @param y y value         
		/// @param label label string         
		////////////////////////////////////////////////////////////////////////////////////////////
		void gplot::label(const double& x, const double& y, const std::string& label)
		{
			std::ostringstream tmp_s;
			tmp_s << "set label \"" << label << "\" at " << x << "," << y;
			send2gp(tmp_s.str());
		}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Set windowtitle.
        /// @param name title string         
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::title(const std::string& name)
        {
            std::ostringstream tmp_s;
            tmp_s << "set title \"" << name << " \" ";
            send2gp(tmp_s.str());
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Set X-axis range.
        /// @param xmin xmin
        /// @param xmax xmax      
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::xlim(const double xmin, const double xmax)
        {
            std::ostringstream tmp_s;
            tmp_s << "set xrange [" << xmin << ":" << xmax << "]";
            send2gp(tmp_s.str());
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Set Y-axis range.
        /// @param ymin ymin
        /// @param ymax ymax  
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::ylim(const double ymin, const double ymax)
        {
            std::ostringstream tmp_s;
            tmp_s << "set yrange [" << ymin << ":" << ymax << "]";
            send2gp(tmp_s.str());
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Plot y vs. x.
        /// @param x x vector
        /// @param y y vector      
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::plot_str2( arma::vec &x, arma::vec &y)
        {
            std::ostringstream tmp_s;
            for(unsigned int n=0;n<y.size();n++)
            {
                tmp_s << x(n) << " " << y(n);
                send2gp(tmp_s.str());
                tmp_s.str(""); // Clear buffer
            }
            send2gp("e");
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Plot y vs. x. with label
        /// @param x     x vector
        /// @param y     y vector      
        /// @param label plot label      
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::plot( arma::vec &x, arma::vec &y, const std::string& label="")
        {
            std::ostringstream tmp_s;
			if (label.empty())
			{ 
				send2gp("set nokey");
			    tmp_s << "plot '-' with " << linestyle;
			}
			else
			{
				send2gp("set key noautotitle");
				tmp_s << "plot '-' title \"" << label << "\" with " << linestyle;
			}            
			send2gp("set grid");
            send2gp(tmp_s.str());
            plot_str2(x,y);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Plot y with label
        /// @param y     y vector      
        /// @param label plot label      
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::plot(arma::vec &y, const std::string& label="")
        {
            arma::vec t(y.size());
            t = arma::linspace(1,y.size(),y.size());
            plot(t,y,label);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Plot y1 and y2 vs. x with label
        /// @param x      x vector
        /// @param y1     y1 vector      
        /// @param y2     y2 vector      
        /// @param label1 plot label 1      
        /// @param label2 plot label 2    
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::plot( arma::vec &x, arma::vec &y1, arma::vec &y2, const std::string& label1, const std::string& label2)
        {
            std::ostringstream tmp_s;
            send2gp("set key noautotitle");
            send2gp("set grid");
            tmp_s << "plot '-' title \"" << label1 << "\" with " << linestyle << ", '-' title \"" << label2 << "\" with " << linestyle;
            send2gp(tmp_s.str());
            plot_str2(x,y1);
            plot_str2(x,y2);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Plot y1 and y2 with label
        /// @param y1     y1 vector      
        /// @param y2     y2 vector      
        /// @param label1 plot label 1      
        /// @param label2 plot label 2    
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::plot(arma::vec &y1, arma::vec &y2, const std::string& label1, const std::string& label2)
        {
            arma::vec t(y1.size());
            t = arma::linspace(1,y1.size(),y1.size());
            plot(t,y1,label1);
            plot(t,y2,label2);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Scatter plot 
        /// @param x    x vector      
        /// @param y    y vector      
        /// @param label plot label      
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::scatter( arma::vec &x, arma::vec &y, const std::string& label="")
        {
            set_linestyle("points");
            plot(x,y,label);
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Scatter plot - dual 
        /// @param x1     x1 vector      
        /// @param y1     y1 vector      
        /// @param label1 plot label1      
        /// @param x2     x2 vector      
        /// @param y2     y2 vector      
        /// @param label2 plot label2      
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::scatter( arma::vec &x1, arma::vec &y1, const std::string& label1,
            arma::vec &x2, arma::vec &y2, const std::string& label2)
        {
            set_linestyle("points");
            std::ostringstream tmp_s;
            send2gp("set key noautotitle");
            send2gp("set grid");
            tmp_s << "plot '-' title \"" << label1 << "\" with " << linestyle << ", '-' title \"" << label2 << "\" with " << linestyle;
            send2gp(tmp_s.str());
            plot_str2(x1,y1);
            plot_str2(x2,y2);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Plot mat as image
        /// @param x     x matrix      
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::image( arma::mat &x)
        {
            std::ostringstream tmp_s;
            xlim(-0.5,x.n_cols-0.5);
            ylim(x.n_rows-0.5,-0.5);
            tmp_s.str(""); // Clear buffer
            tmp_s << "plot '-' matrix with image";
            send2gp(tmp_s.str());
            for(unsigned int r=0;r<x.n_rows;r++)
            {
                tmp_s.str("");  // Clear buffer
                for(unsigned int c=0;c<x.n_cols;c++)
                {
                    tmp_s << x(r,c) << " " ;
                }
                send2gp(tmp_s.str());
            }
            send2gp("e");
            send2gp("e");
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Plot mat as mesh
        /// @param x     x matrix      
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::mesh( arma::mat &x)
        {
            std::ostringstream tmp_s;
            send2gp("unset key");
            send2gp("set hidden3d");
            tmp_s << "splot '-' with lines";
            send2gp(tmp_s.str());

            for(unsigned int r=0;r<x.n_rows;r++)
            {
                for(unsigned int c=0;c<x.n_cols;c++)
                {
                    tmp_s.str("");  // Clear buffer
                    tmp_s << r << " " << c << " "<< x(r,c);
                    send2gp(tmp_s.str());
                }
                send2gp("");
            }
            send2gp("e");
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Plot mat as surf
        /// @param x     x matrix      
        ////////////////////////////////////////////////////////////////////////////////////////////
        void gplot::surf( arma::mat &x)
        {
            send2gp("set pm3d");
            mesh(x);
        }


		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Save plot to file.
		/// @param name filename
		///
		/// Extensions that are supported:
		/// - png
		/// - ps
		/// - eps
		/// - tex
		/// - pdf
		/// - svg
		/// - emf
		/// - gif
		///
		/// \note When 'latex' output is used the '\' must be escaped by '\\\\' e.g set_xlabel("Frequency $\\\\omega = 2 \\\\pi f$")
		////////////////////////////////////////////////////////////////////////////////////////////
		void gplot::set_output(const std::string& name)
		{
			size_t found = name.find_last_of(".");
			std::string ext;
			ext = name.substr(found + 1);
			std::ostringstream tmp_s;

			if (ext.compare("png")==0)
			{
				tmp_s << "set terminal pngcairo enhanced font 'Verdana,10'";
			}
			else if (ext.compare("ps") == 0)
			{
				tmp_s << "set terminal postscript enhanced color";
			}
			else if (ext.compare("eps") == 0)
			{
				tmp_s << "set terminal postscript eps enhanced color";
			}
			else if (ext.compare("tex") == 0)
			{
				tmp_s << "set terminal cairolatex eps color enhanced";
			}
			else if (ext.compare("pdf") == 0)
			{
				tmp_s << "set terminal pdfcairo color enhanced";
			}
			else if (ext.compare("svg") == 0)
			{
				tmp_s << "set terminal svg enhanced";
			}
			else if (ext.compare("emf") == 0)
			{
				tmp_s << "set terminal emf color enhanced";
			}
			else if (ext.compare("gif") == 0)
			{
				tmp_s << "set terminal gif enhanced";
			}
			//else if (ext.compare("jpg") == 0)
			//{
			//	tmp_s << "set terminal jpeg ";
			//}
			else
			{
				tmp_s << "set terminal " << GP_TERM;
			}

			send2gp(tmp_s.str());
			tmp_s.str("");  // Clear buffer
			tmp_s << "set output '" << name << "'";
			send2gp(tmp_s.str());
		}


		////////////////////////////////////////////////////////////////////////////////////////////
		/// \brief Restore output terminal.
		////////////////////////////////////////////////////////////////////////////////////////////
		void gplot::restore_output(void)
		{
			std::ostringstream tmp_s;
			tmp_s << "set terminal " << GP_TERM;
			send2gp(tmp_s.str());
		}
    }; // End Gnuplot Class

  
    /// @}

} // end namespace
#endif
