# ProjCode4YangLiu18Aug

(MatLab) Programs to estimate DCC-MIDAS model with different weighting strategies (I have attached the original paper of the model).
I pushed the code online for backup.



## Update

### March 6, 2019

Hi! I must apology for my mistake: I forgot to cite the original paper in README.
So, I add it here in BibTex:

> @article{COLACITO2011A,
>  title={A component model for dynamic correlations},
>  author={COLACITO and Riccardo and ENGLE and Robert, F and GHYSELS and Eric},
>  journal={Social Science Electronic Publishing},
>  volume={164},
>  number={1},
>  pages={45-59},
>  year={2011},
>}

And, please be careful when using the programs for _Almon_ and _RestrictedBeta_.
Because the theoretical model with these weight stratgies have severe incidental parameter problems (Lancaster, 2000).
The _lessParsUnrestrictedBeta_, which dropped many meaningless parameters, is to solve this problem.
And personally, I recommend _UnrestrictedBeta_. It has much less but enough parameters. And the estimates are stable enough.
(Well, finally I keep _Almon_ & _RestrictedBeta_, after all, this is an entrusted project not my own ideas.)

> Lancaster, T. The incidental parameter problem since 1948. Journal of Econometrics 2000, 95, 391-413.

