//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Unit1.h"
#include <vector>
#include <cmath>
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "VCLTee.TeeSurfa"
#pragma resource "*.dfm"
TForm1 *Form1;
bool stop=true;
    const double _dx = .01, _dy = .01, _dt = (_dx + _dy) / 2. / 100.; //, _dt;

	const size_t N_x = 200;
	const size_t N_y = 100;
	const double _U0 = 1;
	double Re;
	std::vector<std::vector<double> > psi(N_x, std::vector<double>(N_y));
	std::vector<std::vector<double> > e(N_x, std::vector<double>(N_y));
	std::vector<std::vector<double> > u(N_x, std::vector<double>(N_y));
	std::vector<std::vector<double> > v(N_x, std::vector<double>(N_y));
#define BARIER
	// стена, твердая граница
#ifdef BARIER
	const size_t b = static_cast<size_t>(N_y / 2. - N_y * .05);
	const size_t t = static_cast<size_t>(N_y / 2. + N_y * .05);
	const size_t l = static_cast<size_t>(N_x / 5. - N_x * .007);
	const size_t r = static_cast<size_t>(N_x / 5. + N_x * .007);
#endif
    double w = 2. / (1 + sin(M_PI * (_dx + _dy) / 2.));

	size_t iter_MAX = 5000;
	const double eps = 1e-2;

	size_t iter_count = 0;
	size_t iter_best = iter_MAX;
	double w_best = 0;
	double w_step = 15e-3;
	double ct = _dt;
	void start()
	{
      for (size_t i = 0; i < N_x; i++) {
		for (size_t j = 0; j < N_y; j++) {
			psi[i][j] = 0;
			e[i][j] = 0;
			u[i][j] = 0;
			v[i][j] = 0;
		}
	}
    }
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TForm1::StartClick(TObject *Sender)
{
if (stop)
{
   stop = false;
   Start->Caption = "Stop";
}
else
{
   stop = true;
   Start->Caption = "Continue";
}

	while (!stop)
	{
		 ct += _dt;
        { // адвекция диффузия схема против потока
            for (size_t i = 1; i < N_x - 1; i++) {
                for (size_t j = 1; j < N_y - 1; j++) {
#ifdef BARIER
					if ((i >= l && i <= r) && (j >= b && j <= t)) {
                        continue;
                    }
#endif
                    unsigned a = 1, b = 1;
                    if (u[i][j] > 0) {
						a = 0;
                    }
                    if (v[i][j] > 0) {
                        b = 0;
                    }
                    e[i][j] =
						e[i][j] -
                        _dt / _dx * u[i][j] * (1 - a) *
                            (e[i][j] - e[i - 1][j]) -
                        _dt / _dx * a * u[i][j] * (e[i + 1][j] - e[i][j]) -
                        _dt / _dy * v[i][j] * (1 - b) *
                            (e[i][j] - e[i][j - 1]) -
						_dt / _dy * b * v[i][j] * (e[i][j + 1] - e[i][j]) +
                        _dt / Re *
                            ((e[i + 1][j] - 2. * e[i][j] + e[i - 1][j]) / _dx /
                                    _dx +
                                (e[i][j + 1] - 2. * e[i][j] + e[i][j - 1]) /
                                    _dy / _dy);
				}
            }
        }

        double E = 0;
        { //элиптическое
			iter_count = 0;
            double q = pow(_dx / _dy, 2);
            do {
                E = 0;
                for (size_t i = 1; i < N_x - 1; i++) {
                    for (size_t j = 1; j < N_y - 1; j++) {
#ifdef BARIER
                        if ((i >= l && i <= r) && (j >= b && j <= t)) {
                            continue;
                        }
#endif
                        double prev = psi[i][j];
						psi[i][j] = (psi[i + 1][j] + psi[i - 1][j] +
                                        q * (psi[i][j + 1] + psi[i][j - 1])) /
                                        ((2. + 2. * q)) -
                                    _dx * _dx / (2. + 2. * q) * e[i][j];
                        psi[i][j] = prev + w * (psi[i][j] - prev);
                        double er = abs(prev - psi[i][j]);
						if (E < er) {
                            E = er;
                        }
                    }
                }
                iter_count++;
			} while (eps < E && iter_count < iter_MAX);
            if (iter_count <= iter_best) {
                if (iter_count == 1) {
                    iter_best = iter_MAX;
                    w += w_step;
                    if (w > 1.9) {
						w = 1;
                    }
                } else {
                    iter_best = iter_count;
                }
            } else {
				iter_best = iter_MAX;
                w += w_step;
                if (w > 1.9) {
                    w = 1;
                }
            }
		}

        { // пересчет поля скоростей
            for (size_t i = 0; i < N_x; i++) {
                for (size_t j = 1; j < N_y - 1; j++) {
#ifdef BARIER
					if ((i >= l && i <= r) && (j >= b && j <= t)) {
                        continue;
                    }
#endif
                    u[i][j] = (psi[i][j + 1] - psi[i][j]) / _dy;
                }
			}
            for (size_t i = 0; i < N_x - 1; i++) {
                for (size_t j = 1; j < N_y - 1; j++) {
#ifdef BARIER
                    if ((i >= l && i <= r) && (j >= b && j <= t)) {
                        continue;
					}
#endif
                    v[i][j] = -1 * (psi[i + 1][j] - psi[i][j]) / _dx;
                }
            }
        }

#ifdef BARIER
        {
            { // стена, твердая граница
				{ //вихрь на стенках
					for (size_t i = l + 1; i < r; i++) {
						u[i][b] = 0;
						u[i][t] = 0;
						v[i][b] = 0;
						v[i][t] = 0;
						psi[i][b] = 0;
						psi[i][t] = 0;
						e[i][b] = 2. * (psi[i][b - 1]) / _dy / _dy;
						e[i][t] = 2. * (psi[i][t + 1]) / _dy / _dy;
					}
					for (size_t j = b + 1; j < t + 1; j++) {
                        u[l][j] = 0;
                        u[r][j] = 0;
                        v[l][j] = 0;
                        v[r][j] = 0;
						psi[l][j] = 0;
                        psi[r][j] = 0;
                        e[l][j] = 2. * (psi[l - 1][j]) / _dx / _dx;
                        e[r][j] = 2. * (psi[r + 1][j]) / _dx / _dx;
                    }
				}
                { // критические точки на углах
                    e[l][b] = 2. * psi[l][b - 1] / _dy / _dy +
                              2. * psi[l - 1][b] / _dx / _dx;
                    e[r][b] = 2. * psi[r][b - 1] / _dy / _dy +
							  2. * psi[r + 1][b] / _dx / _dx;

                    e[r][t] = 2. * psi[r][t + 1] / _dy / _dy +
                              2. * psi[r + 1][t] / _dx / _dx;
                    e[l][t] = 2. * psi[l][t + 1] / _dy / _dy +
							  2. * psi[l - 1][t] / _dx / _dx;
                }
            }
        }
#endif

        {
            for (size_t i = 0; i < N_x; i++) {
                double U0 = _U0;
                // неограниченая граница верх и низ
				u[i][0] = U0;
                v[i][0] = 0;
                psi[i][0] = psi[i][1] - U0 * _dy;
                e[i][0] = 2. * (psi[i][1] - psi[i][0] - U0 * _dy) / _dy / _dy;

				u[i][N_y - 1] = U0;
                v[i][N_y - 1] = 0;
                psi[i][N_y - 1] = psi[i][N_y - 2] + U0 * _dy;
                e[i][N_y - 1] = 2. *
                                (psi[i][N_y - 2] - psi[i][N_y - 1] + U0 * _dy) /
								_dy / _dy;
            }
        }

        {
			for (size_t j = 1; j < N_y - 1; j++) {
                // градиент функции вихря и тока равен нулю на правой и левой границе
                e[0][j] = e[1][j];
                e[N_x - 1][j] = e[N_x - 2][j];
                psi[0][j] = psi[1][j];
				psi[N_x - 1][j] = psi[N_x - 2][j];
            }
        }

		{ // вывод
				Form1->Series1->Clear();
                for (size_t i = 0; i < N_x; i++) {
                    for (size_t j = 0; j < N_y; j++) {
						Form1->Series1->AddXYZ(i,
							sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j])
							//u[i][j]
							, j);
                    }
                }
                Form1->Caption =
					"Re" + FloatToStr(Re) + " | t:" + FloatToStr(ct);
				Form1->Chart1->Repaint();
            Application->ProcessMessages();
		}
	}
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button1Click(TObject *Sender)
{
   ct=0;
   //Start->Enabled = true;
   Re = StrToFloat(LabeledEdit1->Text);
   start();
   stop = true;
   Start->Caption = "Start";
}
//---------------------------------------------------------------------------

