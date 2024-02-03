      SUBROUTINE DSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE;
      bool               TRYRAC;
      int                IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N;
      double           VL, VU;
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      double             D( * ), E( * ), W( * ), WORK( * );
      double             Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, FOUR, MINRGP;
      const              ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0, MINRGP = 1.0D-3 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LQUERY, VALEIG, WANTZ, ZQUERY, LAESWAP;
      int                I, IBEGIN, IEND, IFIRST, IIL, IINDBL, IINDW, IINDWK, IINFO, IINSPL, IIU, ILAST, IN, INDD, INDE2, INDERR, INDGP, INDGRS, INDWRK, ITMP, ITMP2, J, JBLK, JJ, LIWMIN, LWMIN, NSPLIT, NZCMIN, OFFSET, WBEGIN, WEND;
      double             BIGNUM, CS, EPS, PIVMIN, R1, R2, RMAX, RMIN, RTOL1, RTOL2, SAFMIN, SCALE, SMLNUM, SN, THRESH, TMP, TNRM, WL, WU;
      // ..
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANST;
      // EXTERNAL LSAME, DLAMCH, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLAE2, DLAEV2, DLARRC, DLARRE, DLARRJ, DLARRR, DLARRV, DLASRT, DSCAL, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT


      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      LQUERY = ( ( LWORK.EQ.-1 ).OR.( LIWORK.EQ.-1 ) )
      ZQUERY = ( NZC.EQ.-1 )
      LAESWAP = .FALSE.

      // DSTEMR needs WORK of size 6*N, IWORK of size 3*N.
      // In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N.
      // Furthermore, DLARRV needs WORK of size 12*N, IWORK of size 7*N.
      if ( WANTZ ) {
         LWMIN = 18*N
         LIWMIN = 10*N
      } else {
         // need less workspace if only the eigenvalues are wanted
         LWMIN = 12*N
         LIWMIN = 8*N
      ENDIF

      WL = ZERO
      WU = ZERO
      IIL = 0
      IIU = 0
      NSPLIT = 0

      if ( VALEIG ) {
         // We do not reference VL, VU in the cases RANGE = 'I','A'
         // The interval (WL, WU] contains all the wanted eigenvalues.
         // It is either given by the user or computed in DLARRE.
         WL = VL
         WU = VU
      ELSEIF( INDEIG ) THEN
         // We do not reference IL, IU in the cases RANGE = 'V','A'
         IIL = IL
         IIU = IU
      ENDIF

      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( VALEIG .AND. N.GT.0 .AND. WU.LE.WL ) {
         INFO = -7
      } else if ( INDEIG .AND. ( IIL.LT.1 .OR. IIL.GT.N ) ) {
         INFO = -8
      } else if ( INDEIG .AND. ( IIU.LT.IIL .OR. IIU.GT.N ) ) {
         INFO = -9
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -13
      } else if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
         INFO = -17
      } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
         INFO = -19
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )

      if ( INFO.EQ.0 ) {
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN

         if ( WANTZ .AND. ALLEIG ) {
            NZCMIN = N
         } else if ( WANTZ .AND. VALEIG ) {
            dlarrc('T', N, VL, VU, D, E, SAFMIN, NZCMIN, ITMP, ITMP2, INFO );
         } else if ( WANTZ .AND. INDEIG ) {
            NZCMIN = IIU-IIL+1
         } else {
            // WANTZ .EQ. FALSE.
            NZCMIN = 0
         ENDIF
         if ( ZQUERY .AND. INFO.EQ.0 ) {
            Z( 1,1 ) = NZCMIN
         } else if ( NZC.LT.NZCMIN .AND. .NOT.ZQUERY ) {
            INFO = -14
         }
      }

      if ( INFO.NE.0 ) {

         xerbla('DSTEMR', -INFO );

         RETURN
      } else if ( LQUERY .OR. ZQUERY ) {
         RETURN
      }

      // Handle N = 0, 1, and 2 cases immediately

      M = 0
      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         if ( ALLEIG .OR. INDEIG ) {
            M = 1
            W( 1 ) = D( 1 )
         } else {
            if ( WL.LT.D( 1 ) .AND. WU.GE.D( 1 ) ) {
               M = 1
               W( 1 ) = D( 1 )
            }
         }
         if ( WANTZ.AND.(.NOT.ZQUERY) ) {
            Z( 1, 1 ) = ONE
            ISUPPZ(1) = 1
            ISUPPZ(2) = 1
         }
         RETURN
      }

      if ( N.EQ.2 ) {
         if ( .NOT.WANTZ ) {
            dlae2(D(1), E(1), D(2), R1, R2 );
         } else if ( WANTZ.AND.(.NOT.ZQUERY) ) {
            dlaev2(D(1), E(1), D(2), R1, R2, CS, SN );
         }
         // D/S/LAE2 and D/S/LAEV2 outputs satisfy |R1| >= |R2|. However,
         // the following code requires R1 >= R2. Hence, we correct
         // the order of R1, R2, CS, SN if R1 < R2 before further processing.
         if ( R1.LT.R2 ) {
            E(2) = R1
            R1 = R2
            R2 = E(2)
            LAESWAP = .TRUE.
         ENDIF
         if ( ALLEIG.OR. (VALEIG.AND.(R2.GT.WL).AND. (R2.LE.WU)).OR. (INDEIG.AND.(IIL.EQ.1)) ) {
            M = M+1
            W( M ) = R2
            if ( WANTZ.AND.(.NOT.ZQUERY) ) {
               if ( LAESWAP ) {
                  Z( 1, M ) = CS
                  Z( 2, M ) = SN
               } else {
                  Z( 1, M ) = -SN
                  Z( 2, M ) = CS
               ENDIF
               // Note: At most one of SN and CS can be zero.
               if (SN.NE.ZERO) {
                  if (CS.NE.ZERO) {
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 2
                  } else {
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 1
                  }
               } else {
                  ISUPPZ(2*M-1) = 2
                  ISUPPZ(2*M) = 2
               }
            ENDIF
         ENDIF
         if ( ALLEIG.OR. (VALEIG.AND.(R1.GT.WL).AND. (R1.LE.WU)).OR. (INDEIG.AND.(IIU.EQ.2)) ) {
            M = M+1
            W( M ) = R1
            if ( WANTZ.AND.(.NOT.ZQUERY) ) {
               if ( LAESWAP ) {
                  Z( 1, M ) = -SN
                  Z( 2, M ) = CS
               } else {
                  Z( 1, M ) = CS
                  Z( 2, M ) = SN
               ENDIF
               // Note: At most one of SN and CS can be zero.
               if (SN.NE.ZERO) {
                  if (CS.NE.ZERO) {
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 2
                  } else {
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 1
                  }
               } else {
                  ISUPPZ(2*M-1) = 2
                  ISUPPZ(2*M) = 2
               }
            ENDIF
         ENDIF

      } else {

      // Continue with general N

         INDGRS = 1
         INDERR = 2*N + 1
         INDGP = 3*N + 1
         INDD = 4*N + 1
         INDE2 = 5*N + 1
         INDWRK = 6*N + 1

         IINSPL = 1
         IINDBL = N + 1
         IINDW = 2*N + 1
         IINDWK = 3*N + 1

         // Scale matrix to allowable range, if necessary.
         // The allowable range is related to the PIVMIN parameter; see the
         // comments in DLARRD.  The preference for scaling small values
         // up is heuristic; we expect users' matrices not to be close to the
         // RMAX threshold.

         SCALE = ONE
         TNRM = DLANST( 'M', N, D, E )
         if ( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) {
            SCALE = RMIN / TNRM
         } else if ( TNRM.GT.RMAX ) {
            SCALE = RMAX / TNRM
         }
         if ( SCALE.NE.ONE ) {
            dscal(N, SCALE, D, 1 );
            dscal(N-1, SCALE, E, 1 );
            TNRM = TNRM*SCALE
            if ( VALEIG ) {
               // If eigenvalues in interval have to be found,
               // scale (WL, WU] accordingly
               WL = WL*SCALE
               WU = WU*SCALE
            ENDIF
         }

         // Compute the desired eigenvalues of the tridiagonal after splitting
         // into smaller subblocks if the corresponding off-diagonal elements
         // are small
         // THRESH is the splitting parameter for DLARRE
         // A negative THRESH forces the old splitting criterion based on the
         // size of the off-diagonal. A positive THRESH switches to splitting
         // which preserves relative accuracy.

         if ( TRYRAC ) {
            // Test whether the matrix warrants the more expensive relative approach.
            dlarrr(N, D, E, IINFO );
         } else {
            // The user does not care about relative accurately eigenvalues
            IINFO = -1
         ENDIF
         // Set the splitting criterion
         if (IINFO.EQ.0) {
            THRESH = EPS
         } else {
            THRESH = -EPS
            // relative accuracy is desired but T does not guarantee it
            TRYRAC = .FALSE.
         ENDIF

         if ( TRYRAC ) {
            // Copy original diagonal, needed to guarantee relative accuracy
            dcopy(N,D,1,WORK(INDD),1);
         ENDIF
         // Store the squares of the offdiagonal values of T
         DO 5 J = 1, N-1
            WORK( INDE2+J-1 ) = E(J)**2
         } // 5

         // Set the tolerance parameters for bisection
         if ( .NOT.WANTZ ) {
            // DLARRE computes the eigenvalues to full precision.
            RTOL1 = FOUR * EPS
            RTOL2 = FOUR * EPS
         } else {
            // DLARRE computes the eigenvalues to less than full precision.
            // DLARRV will refine the eigenvalue approximations, and we can
            // need less accurate initial bisection in DLARRE.
            // Note: these settings do only affect the subset case and DLARRE
            RTOL1 = SQRT(EPS)
            RTOL2 = MAX( SQRT(EPS)*5.0D-3, FOUR * EPS )
         ENDIF
         dlarre(RANGE, N, WL, WU, IIL, IIU, D, E, WORK(INDE2), RTOL1, RTOL2, THRESH, NSPLIT, IWORK( IINSPL ), M, W, WORK( INDERR ), WORK( INDGP ), IWORK( IINDBL ), IWORK( IINDW ), WORK( INDGRS ), PIVMIN, WORK( INDWRK ), IWORK( IINDWK ), IINFO );
         if ( IINFO.NE.0 ) {
            INFO = 10 + ABS( IINFO )
            RETURN
         }
         // Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired
         // part of the spectrum. All desired eigenvalues are contained in
         // (WL,WU]


         if ( WANTZ ) {

            // Compute the desired eigenvectors corresponding to the computed
            // eigenvalues

            dlarrv(N, WL, WU, D, E, PIVMIN, IWORK( IINSPL ), M, 1, M, MINRGP, RTOL1, RTOL2, W, WORK( INDERR ), WORK( INDGP ), IWORK( IINDBL ), IWORK( IINDW ), WORK( INDGRS ), Z, LDZ, ISUPPZ, WORK( INDWRK ), IWORK( IINDWK ), IINFO );
            if ( IINFO.NE.0 ) {
               INFO = 20 + ABS( IINFO )
               RETURN
            }
         } else {
            // DLARRE computes eigenvalues of the (shifted) root representation
            // DLARRV returns the eigenvalues of the unshifted matrix.
            // However, if the eigenvectors are not desired by the user, we need
            // to apply the corresponding shifts from DLARRE to obtain the
            // eigenvalues of the original matrix.
            for (J = 1; J <= M; J++) { // 20
               ITMP = IWORK( IINDBL+J-1 )
               W( J ) = W( J ) + E( IWORK( IINSPL+ITMP-1 ) )
            } // 20
         }


         if ( TRYRAC ) {
            // Refine computed eigenvalues so that they are relatively accurate
            // with respect to the original matrix T.
            IBEGIN = 1
            WBEGIN = 1
            DO 39  JBLK = 1, IWORK( IINDBL+M-1 )
               IEND = IWORK( IINSPL+JBLK-1 )
               IN = IEND - IBEGIN + 1
               WEND = WBEGIN - 1
               // check if any eigenvalues have to be refined in this block
               } // 36
               if ( WEND.LT.M ) {
                  if ( IWORK( IINDBL+WEND ).EQ.JBLK ) {
                     WEND = WEND + 1
                     GO TO 36
                  }
               }
               if ( WEND.LT.WBEGIN ) {
                  IBEGIN = IEND + 1
                  GO TO 39
               }

               OFFSET = IWORK(IINDW+WBEGIN-1)-1
               IFIRST = IWORK(IINDW+WBEGIN-1)
               ILAST = IWORK(IINDW+WEND-1)
               RTOL2 = FOUR * EPS
               dlarrj(IN, WORK(INDD+IBEGIN-1), WORK(INDE2+IBEGIN-1), IFIRST, ILAST, RTOL2, OFFSET, W(WBEGIN), WORK( INDERR+WBEGIN-1 ), WORK( INDWRK ), IWORK( IINDWK ), PIVMIN, TNRM, IINFO );
               IBEGIN = IEND + 1
               WBEGIN = WEND + 1
            } // 39
         ENDIF

         // If matrix was scaled, then rescale eigenvalues appropriately.

         if ( SCALE.NE.ONE ) {
            dscal(M, ONE / SCALE, W, 1 );
         }

      }


      // If eigenvalues are not in increasing order, then sort them,
      // possibly along with eigenvectors.

      if ( NSPLIT.GT.1 .OR. N.EQ.2 ) {
         if ( .NOT. WANTZ ) {
            dlasrt('I', M, W, IINFO );
            if ( IINFO.NE.0 ) {
               INFO = 3
               RETURN
            }
         } else {
            DO 60 J = 1, M - 1
               I = 0
               TMP = W( J )
               DO 50 JJ = J + 1, M
                  if ( W( JJ ).LT.TMP ) {
                     I = JJ
                     TMP = W( JJ )
                  }
               } // 50
               if ( I.NE.0 ) {
                  W( I ) = W( J )
                  W( J ) = TMP
                  if ( WANTZ ) {
                     dswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
                     ITMP = ISUPPZ( 2*I-1 )
                     ISUPPZ( 2*I-1 ) = ISUPPZ( 2*J-1 )
                     ISUPPZ( 2*J-1 ) = ITMP
                     ITMP = ISUPPZ( 2*I )
                     ISUPPZ( 2*I ) = ISUPPZ( 2*J )
                     ISUPPZ( 2*J ) = ITMP
                  }
               }
            } // 60
         }
      ENDIF


      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of DSTEMR

      }
