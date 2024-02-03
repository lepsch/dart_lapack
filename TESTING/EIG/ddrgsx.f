      void ddrgsx(NSIZE, NCMAX, THRESH, NIN, NOUT, A, LDA, B, AI, BI, Z, Q, ALPHAR, ALPHAI, BETA, C, LDC, S, WORK, LWORK, IWORK, LIWORK, BWORK, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDC, LIWORK, LWORK, NCMAX, NIN, NOUT, NSIZE;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      double             A( LDA, * ), AI( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDA, * ), BETA( * ), BI( LDA, * ), C( LDC, * ), Q( LDA, * ), S( * ), WORK( * ), Z( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0, ONE = 1.0, TEN = 1.0e+1 ;
      // ..
      // .. Local Scalars ..
      bool               ILABAD;
      String             SENSE;
      int                BDSPAC, I, I1, IFUNC, IINFO, J, LINFO, MAXWRK, MINWRK, MM, MN2, NERRS, NPTKNT, NTEST, NTESTT, PRTYPE, QBA, QBB;
      double             ABNRM, BIGNUM, DIFTRU, PLTRU, SMLNUM, TEMP1, TEMP2, THRSH2, ULP, ULPINV, WEIGHT;
      // ..
      // .. Local Arrays ..
      double             DIFEST( 2 ), PL( 2 ), RESULT( 10 );
      // ..
      // .. External Functions ..
      //- bool               DLCTSX;
      //- int                ILAENV;
      //- double             DLAMCH, DLANGE;
      // EXTERNAL DLCTSX, ILAENV, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, DGESVD, DGET51, DGET53, DGGESX, DLACPY, DLAKF2, DLASET, DLATM5, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Scalars in Common ..
      bool               FS;
      int                K, M, MPLUSN, N;
      // ..
      // .. Common blocks ..
      // COMMON / MN / M, N, MPLUSN, K, FS
      // ..
      // .. Executable Statements ..

      // Check for errors

      if ( NSIZE < 0 ) {
         INFO = -1;
      } else if ( THRESH < ZERO ) {
         INFO = -2;
      } else if ( NIN <= 0 ) {
         INFO = -3;
      } else if ( NOUT <= 0 ) {
         INFO = -4;
      } else if ( LDA < 1 || LDA < NSIZE ) {
         INFO = -6;
      } else if ( LDC < 1 || LDC < NSIZE*NSIZE / 2 ) {
         INFO = -17;
      } else if ( LIWORK < NSIZE+6 ) {
         INFO = -21;
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      MINWRK = 1;
      if ( INFO == 0 && LWORK >= 1 ) {
         MINWRK = max( 10*( NSIZE+1 ), 5*NSIZE*NSIZE / 2 );

         // workspace for sggesx

         MAXWRK = 9*( NSIZE+1 ) + NSIZE* ILAENV( 1, 'DGEQRF', ' ', NSIZE, 1, NSIZE, 0 )          MAXWRK = max( MAXWRK, 9*( NSIZE+1 )+NSIZE* ILAENV( 1, 'DORGQR', ' ', NSIZE, 1, NSIZE, -1 ) );

         // workspace for dgesvd

         BDSPAC = 5*NSIZE*NSIZE / 2;
         MAXWRK = max( MAXWRK, 3*NSIZE*NSIZE / 2+NSIZE*NSIZE* ILAENV( 1, 'DGEBRD', ' ', NSIZE*NSIZE / 2, NSIZE*NSIZE / 2, -1, -1 ) );
         MAXWRK = max( MAXWRK, BDSPAC );

         MAXWRK = max( MAXWRK, MINWRK );

         WORK( 1 ) = MAXWRK;
      }

      if (LWORK < MINWRK) INFO = -19;

      if ( INFO != 0 ) {
         xerbla('DDRGSX', -INFO );
         return;
      }

      // Important constants

      ULP = DLAMCH( 'P' );
      ULPINV = ONE / ULP;
      SMLNUM = DLAMCH( 'S' ) / ULP;
      BIGNUM = ONE / SMLNUM;
      THRSH2 = TEN*THRESH;
      NTESTT = 0;
      NERRS = 0;

      // Go to the tests for read-in matrix pairs

      IFUNC = 0;
      if (NSIZE == 0) GO TO 70;

      // Test the built-in matrix pairs.
      // Loop over different functions (IFUNC) of DGGESX, types (PRTYPE)
      // of test matrices, different size (M+N)

      PRTYPE = 0;
      QBA = 3;
      QBB = 4;
      WEIGHT = sqrt( ULP );

      for (IFUNC = 0; IFUNC <= 3; IFUNC++) { // 60
         for (PRTYPE = 1; PRTYPE <= 5; PRTYPE++) { // 50
            for (M = 1; M <= NSIZE - 1; M++) { // 40
               for (N = 1; N <= NSIZE - M; N++) { // 30

                  WEIGHT = ONE / WEIGHT;
                  MPLUSN = M + N;

                  // Generate test matrices

                  FS = true;
                  K = 0;

                  dlaset('Full', MPLUSN, MPLUSN, ZERO, ZERO, AI, LDA );
                  dlaset('Full', MPLUSN, MPLUSN, ZERO, ZERO, BI, LDA );

                  dlatm5(PRTYPE, M, N, AI, LDA, AI( M+1, M+1 ), LDA, AI( 1, M+1 ), LDA, BI, LDA, BI( M+1, M+1 ), LDA, BI( 1, M+1 ), LDA, Q, LDA, Z, LDA, WEIGHT, QBA, QBB );

                  // Compute the Schur factorization and swapping the
                  // m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.
                  // Swapping is accomplished via the function DLCTSX
                  // which is supplied below.

                  if ( IFUNC == 0 ) {
                     SENSE = 'N';
                  } else if ( IFUNC == 1 ) {
                     SENSE = 'E';
                  } else if ( IFUNC == 2 ) {
                     SENSE = 'V';
                  } else if ( IFUNC == 3 ) {
                     SENSE = 'B';
                  }

                  dlacpy('Full', MPLUSN, MPLUSN, AI, LDA, A, LDA );
                  dlacpy('Full', MPLUSN, MPLUSN, BI, LDA, B, LDA );

                  dggesx('V', 'V', 'S', DLCTSX, SENSE, MPLUSN, AI, LDA, BI, LDA, MM, ALPHAR, ALPHAI, BETA, Q, LDA, Z, LDA, PL, DIFEST, WORK, LWORK, IWORK, LIWORK, BWORK, LINFO );

                  if ( LINFO != 0 && LINFO != MPLUSN+2 ) {
                     RESULT( 1 ) = ULPINV;
                     WRITE( NOUT, FMT = 9999 )'DGGESX', LINFO, MPLUSN, PRTYPE;
                     INFO = LINFO;
                     GO TO 30;
                  }

                  // Compute the norm(A, B)

                  dlacpy('Full', MPLUSN, MPLUSN, AI, LDA, WORK, MPLUSN );
                  dlacpy('Full', MPLUSN, MPLUSN, BI, LDA, WORK( MPLUSN*MPLUSN+1 ), MPLUSN )                   ABNRM = DLANGE( 'Fro', MPLUSN, 2*MPLUSN, WORK, MPLUSN, WORK );

                  // Do tests (1) to (4)

                  dget51(1, MPLUSN, A, LDA, AI, LDA, Q, LDA, Z, LDA, WORK, RESULT( 1 ) );
                  dget51(1, MPLUSN, B, LDA, BI, LDA, Q, LDA, Z, LDA, WORK, RESULT( 2 ) );
                  dget51(3, MPLUSN, B, LDA, BI, LDA, Q, LDA, Q, LDA, WORK, RESULT( 3 ) );
                  dget51(3, MPLUSN, B, LDA, BI, LDA, Z, LDA, Z, LDA, WORK, RESULT( 4 ) );
                  NTEST = 4;

                  // Do tests (5) and (6): check Schur form of A and
                  // compare eigenvalues with diagonals.

                  TEMP1 = ZERO;
                  RESULT( 5 ) = ZERO;
                  RESULT( 6 ) = ZERO;

                  for (J = 1; J <= MPLUSN; J++) { // 10
                     ILABAD = false;
                     if ( ALPHAI( J ) == ZERO ) {
                        TEMP2 = ( ABS( ALPHAR( J )-AI( J, J ) ) / max( SMLNUM, ( ALPHAR( J ) ).abs(), ( AI( J, J ) ) ).abs()+ ABS( BETA( J )-BI( J, J ) ) / max( SMLNUM, ( BETA( J ) ).abs(), ( BI( J, J ) ) ) ).abs() / ULP;
                        if ( J < MPLUSN ) {
                           if ( AI( J+1, J ) != ZERO ) {
                              ILABAD = true;
                              RESULT( 5 ) = ULPINV;
                           }
                        }
                        if ( J > 1 ) {
                           if ( AI( J, J-1 ) != ZERO ) {
                              ILABAD = true;
                              RESULT( 5 ) = ULPINV;
                           }
                        }
                     } else {
                        if ( ALPHAI( J ) > ZERO ) {
                           I1 = J;
                        } else {
                           I1 = J - 1;
                        }
                        if ( I1 <= 0 || I1 >= MPLUSN ) {
                           ILABAD = true;
                        } else if ( I1 < MPLUSN-1 ) {
                           if ( AI( I1+2, I1+1 ) != ZERO ) {
                              ILABAD = true;
                              RESULT( 5 ) = ULPINV;
                           }
                        } else if ( I1 > 1 ) {
                           if ( AI( I1, I1-1 ) != ZERO ) {
                              ILABAD = true;
                              RESULT( 5 ) = ULPINV;
                           }
                        }
                        if ( !ILABAD ) {
                           dget53(AI( I1, I1 ), LDA, BI( I1, I1 ), LDA, BETA( J ), ALPHAR( J ), ALPHAI( J ), TEMP2, IINFO );
                           if ( IINFO >= 3 ) {
                              WRITE( NOUT, FMT = 9997 )IINFO, J, MPLUSN, PRTYPE;
                              INFO = ( IINFO ).abs();
                           }
                        } else {
                           TEMP2 = ULPINV;
                        }
                     }
                     TEMP1 = max( TEMP1, TEMP2 );
                     if ( ILABAD ) {
                        WRITE( NOUT, FMT = 9996 )J, MPLUSN, PRTYPE;
                     }
                  } // 10
                  RESULT( 6 ) = TEMP1;
                  NTEST = NTEST + 2;

                  // Test (7) (if sorting worked)

                  RESULT( 7 ) = ZERO;
                  if ( LINFO == MPLUSN+3 ) {
                     RESULT( 7 ) = ULPINV;
                  } else if ( MM != N ) {
                     RESULT( 7 ) = ULPINV;
                  }
                  NTEST = NTEST + 1;

                  // Test (8): compare the estimated value DIF and its
                  // value. first, compute the exact DIF.

                  RESULT( 8 ) = ZERO;
                  MN2 = MM*( MPLUSN-MM )*2;
                  if ( IFUNC >= 2 && MN2 <= NCMAX*NCMAX ) {

                     // Note: for either following two causes, there are
                     // almost same number of test cases fail the test.

                     dlakf2(MM, MPLUSN-MM, AI, LDA, AI( MM+1, MM+1 ), BI, BI( MM+1, MM+1 ), C, LDC );

                     dgesvd('N', 'N', MN2, MN2, C, LDC, S, WORK, 1, WORK( 2 ), 1, WORK( 3 ), LWORK-2, INFO );
                     DIFTRU = S( MN2 );

                     if ( DIFEST( 2 ) == ZERO ) {
                        if (DIFTRU > ABNRM*ULP) RESULT( 8 ) = ULPINV;
                     } else if ( DIFTRU == ZERO ) {
                        if( DIFEST( 2 ) > ABNRM*ULP ) RESULT( 8 ) = ULPINV;
                     ELSE IF( ( DIFTRU > THRSH2*DIFEST( 2 ) ) || ( DIFTRU*THRSH2 < DIFEST( 2 ) ) ) {
                        RESULT( 8 ) = max( DIFTRU / DIFEST( 2 ), DIFEST( 2 ) / DIFTRU );
                     }
                     NTEST = NTEST + 1;
                  }

                  // Test (9)

                  RESULT( 9 ) = ZERO;
                  if ( LINFO == ( MPLUSN+2 ) ) {
                     if (DIFTRU > ABNRM*ULP) RESULT( 9 ) = ULPINV;
                     if( ( IFUNC > 1 ) && ( DIFEST( 2 ) != ZERO ) ) RESULT( 9 ) = ULPINV;
                     IF( ( IFUNC == 1 ) && ( PL( 1 ) != ZERO ) ) RESULT( 9 ) = ULPINV;
                     NTEST = NTEST + 1;
                  }

                  NTESTT = NTESTT + NTEST;

                  // Print out tests which fail.

                  for (J = 1; J <= 9; J++) { // 20
                     if ( RESULT( J ) >= THRESH ) {

                        // If this is the first test to fail,
                        // print a header to the data file.

                        if ( NERRS == 0 ) {
                           WRITE( NOUT, FMT = 9995 )'DGX';

                           // Matrix types

                           WRITE( NOUT, FMT = 9993 );

                           // Tests performed

                           WRITE( NOUT, FMT = 9992 )'orthogonal', '''', 'transpose', ( '''', I = 1, 4 );

                        }
                        NERRS = NERRS + 1;
                        if ( RESULT( J ) < 10000.0 ) {
                           WRITE( NOUT, FMT = 9991 )MPLUSN, PRTYPE, WEIGHT, M, J, RESULT( J );
                        } else {
                           WRITE( NOUT, FMT = 9990 )MPLUSN, PRTYPE, WEIGHT, M, J, RESULT( J );
                        }
                     }
                  } // 20

               } // 30
            } // 40
         } // 50
      } // 60

      GO TO 150;

      } // 70

      // Read in data from file to check accuracy of condition estimation
      // Read input data until N=0

      NPTKNT = 0;

      } // 80
      READ( NIN, FMT = *, END = 140 )MPLUSN;
      if (MPLUSN == 0) GO TO 140;
      READ( NIN, FMT = *, END = 140 )N;
      for (I = 1; I <= MPLUSN; I++) { // 90
         READ( NIN, FMT = * )( AI( I, J ), J = 1, MPLUSN );
      } // 90
      for (I = 1; I <= MPLUSN; I++) { // 100
         READ( NIN, FMT = * )( BI( I, J ), J = 1, MPLUSN );
      } // 100
      READ( NIN, FMT = * )PLTRU, DIFTRU;

      NPTKNT = NPTKNT + 1;
      FS = true;
      K = 0;
      M = MPLUSN - N;

      dlacpy('Full', MPLUSN, MPLUSN, AI, LDA, A, LDA );
      dlacpy('Full', MPLUSN, MPLUSN, BI, LDA, B, LDA );

      // Compute the Schur factorization while swapping the
      // m-by-m (1,1)-blocks with n-by-n (2,2)-blocks.

      dggesx('V', 'V', 'S', DLCTSX, 'B', MPLUSN, AI, LDA, BI, LDA, MM, ALPHAR, ALPHAI, BETA, Q, LDA, Z, LDA, PL, DIFEST, WORK, LWORK, IWORK, LIWORK, BWORK, LINFO );

      if ( LINFO != 0 && LINFO != MPLUSN+2 ) {
         RESULT( 1 ) = ULPINV;
         WRITE( NOUT, FMT = 9998 )'DGGESX', LINFO, MPLUSN, NPTKNT;
         GO TO 130;
      }

      // Compute the norm(A, B)
         // (should this be norm of (A,B) or (AI,BI)?)

      dlacpy('Full', MPLUSN, MPLUSN, AI, LDA, WORK, MPLUSN );
      dlacpy('Full', MPLUSN, MPLUSN, BI, LDA, WORK( MPLUSN*MPLUSN+1 ), MPLUSN );
      ABNRM = DLANGE( 'Fro', MPLUSN, 2*MPLUSN, WORK, MPLUSN, WORK );

      // Do tests (1) to (4)

      dget51(1, MPLUSN, A, LDA, AI, LDA, Q, LDA, Z, LDA, WORK, RESULT( 1 ) )       CALL DGET51( 1, MPLUSN, B, LDA, BI, LDA, Q, LDA, Z, LDA, WORK, RESULT( 2 ) )       CALL DGET51( 3, MPLUSN, B, LDA, BI, LDA, Q, LDA, Q, LDA, WORK, RESULT( 3 ) )       CALL DGET51( 3, MPLUSN, B, LDA, BI, LDA, Z, LDA, Z, LDA, WORK, RESULT( 4 ) );

      // Do tests (5) and (6): check Schur form of A and compare
      // eigenvalues with diagonals.

      NTEST = 6;
      TEMP1 = ZERO;
      RESULT( 5 ) = ZERO;
      RESULT( 6 ) = ZERO;

      for (J = 1; J <= MPLUSN; J++) { // 110
         ILABAD = false;
         if ( ALPHAI( J ) == ZERO ) {
            TEMP2 = ( ABS( ALPHAR( J )-AI( J, J ) ) / max( SMLNUM, ( ALPHAR( J ) ).abs(), ( AI( J, J ) ) ).abs()+ABS( BETA( J )-BI( J, J ) ) / max( SMLNUM, ( BETA( J ) ).abs(), ( BI( J, J ) ) ) ).abs() / ULP;
            if ( J < MPLUSN ) {
               if ( AI( J+1, J ) != ZERO ) {
                  ILABAD = true;
                  RESULT( 5 ) = ULPINV;
               }
            }
            if ( J > 1 ) {
               if ( AI( J, J-1 ) != ZERO ) {
                  ILABAD = true;
                  RESULT( 5 ) = ULPINV;
               }
            }
         } else {
            if ( ALPHAI( J ) > ZERO ) {
               I1 = J;
            } else {
               I1 = J - 1;
            }
            if ( I1 <= 0 || I1 >= MPLUSN ) {
               ILABAD = true;
            } else if ( I1 < MPLUSN-1 ) {
               if ( AI( I1+2, I1+1 ) != ZERO ) {
                  ILABAD = true;
                  RESULT( 5 ) = ULPINV;
               }
            } else if ( I1 > 1 ) {
               if ( AI( I1, I1-1 ) != ZERO ) {
                  ILABAD = true;
                  RESULT( 5 ) = ULPINV;
               }
            }
            if ( !ILABAD ) {
               dget53(AI( I1, I1 ), LDA, BI( I1, I1 ), LDA, BETA( J ), ALPHAR( J ), ALPHAI( J ), TEMP2, IINFO );
               if ( IINFO >= 3 ) {
                  WRITE( NOUT, FMT = 9997 )IINFO, J, MPLUSN, NPTKNT;
                  INFO = ( IINFO ).abs();
               }
            } else {
               TEMP2 = ULPINV;
            }
         }
         TEMP1 = max( TEMP1, TEMP2 );
         if ( ILABAD ) {
            WRITE( NOUT, FMT = 9996 )J, MPLUSN, NPTKNT;
         }
      } // 110
      RESULT( 6 ) = TEMP1;

      // Test (7) (if sorting worked)  <--------- need to be checked.

      NTEST = 7;
      RESULT( 7 ) = ZERO;
      if (LINFO == MPLUSN+3) RESULT( 7 ) = ULPINV;

      // Test (8): compare the estimated value of DIF and its true value.

      NTEST = 8;
      RESULT( 8 ) = ZERO;
      if ( DIFEST( 2 ) == ZERO ) {
         if (DIFTRU > ABNRM*ULP) RESULT( 8 ) = ULPINV;
      } else if ( DIFTRU == ZERO ) {
         if ( DIFEST( 2 ) > ABNRM*ULP ) RESULT( 8 ) = ULPINV;
      ELSE IF( ( DIFTRU > THRSH2*DIFEST( 2 ) ) || ( DIFTRU*THRSH2 < DIFEST( 2 ) ) ) {
         RESULT( 8 ) = max( DIFTRU / DIFEST( 2 ), DIFEST( 2 ) / DIFTRU );
      }

      // Test (9)

      NTEST = 9;
      RESULT( 9 ) = ZERO;
      if ( LINFO == ( MPLUSN+2 ) ) {
         if (DIFTRU > ABNRM*ULP) RESULT( 9 ) = ULPINV;
         if( ( IFUNC > 1 ) && ( DIFEST( 2 ) != ZERO ) ) RESULT( 9 ) = ULPINV;
         IF( ( IFUNC == 1 ) && ( PL( 1 ) != ZERO ) ) RESULT( 9 ) = ULPINV;
      }

      // Test (10): compare the estimated value of PL and it true value.

      NTEST = 10;
      RESULT( 10 ) = ZERO;
      if ( PL( 1 ) == ZERO ) {
         if (PLTRU > ABNRM*ULP) RESULT( 10 ) = ULPINV;
      } else if ( PLTRU == ZERO ) {
         if ( PL( 1 ) > ABNRM*ULP ) RESULT( 10 ) = ULPINV;
      ELSE IF( ( PLTRU > THRESH*PL( 1 ) ) || ( PLTRU*THRESH < PL( 1 ) ) ) {
         RESULT( 10 ) = ULPINV;
      }

      NTESTT = NTESTT + NTEST;

      // Print out tests which fail.

      for (J = 1; J <= NTEST; J++) { // 120
         if ( RESULT( J ) >= THRESH ) {

            // If this is the first test to fail,
            // print a header to the data file.

            if ( NERRS == 0 ) {
               WRITE( NOUT, FMT = 9995 )'DGX';

               // Matrix types

               WRITE( NOUT, FMT = 9994 );

               // Tests performed

               WRITE( NOUT, FMT = 9992 )'orthogonal', '''', 'transpose', ( '''', I = 1, 4 );

            }
            NERRS = NERRS + 1;
            if ( RESULT( J ) < 10000.0 ) {
               WRITE( NOUT, FMT = 9989 )NPTKNT, MPLUSN, J, RESULT( J );
            } else {
               WRITE( NOUT, FMT = 9988 )NPTKNT, MPLUSN, J, RESULT( J );
            }
         }

      } // 120

      } // 130
      GO TO 80;
      } // 140

      } // 150

      // Summary

      alasvm('DGX', NOUT, NERRS, NTESTT, 0 );

      WORK( 1 ) = MAXWRK;

      return;

 9999 FORMAT( ' DDRGSX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ')' );

 9998 FORMAT( ' DDRGSX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', Input Example #', I2, ')' );

 9997 FORMAT( ' DDRGSX: DGET53 returned INFO=', I1, ' for eigenvalue ', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ')' );

 9996 FORMAT( ' DDRGSX: S not in Schur form at eigenvalue ', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ')' );

 9995 FORMAT( / 1X, A3, ' -- Real Expert Generalized Schur form', ' problem driver' );

 9994 FORMAT( 'Input Example' );

 9993 FORMAT( ' Matrix types: ', / '  1:  A is a block diagonal matrix of Jordan blocks ', 'and B is the identity ', / '      matrix, ', / '  2:  A and B are upper triangular matrices, ', / '  3:  A and B are as type 2, but each second diagonal ', 'block in A_11 and ', / '      each third diagonal block in A_22 are 2x2 blocks,', / '  4:  A and B are block diagonal matrices, ', / '  5:  (A,B) has potentially close or common ', 'eigenvalues.', / );

 9992 FORMAT( / ' Tests performed:  (S is Schur, T is triangular, ', 'Q and Z are ', A, ',', / 19X, ' a is alpha, b is beta, and ', A, ' means ', A, '.)', / '  1 = | A - Q S Z', A, ' | / ( |A| n ulp )      2 = | B - Q T Z', A, ' | / ( |B| n ulp )', / '  3 = | I - QQ', A, ' | / ( n ulp )             4 = | I - ZZ', A, ' | / ( n ulp )', / '  5 = 1/ULP  if A is not in ', 'Schur form S', / '  6 = difference between (alpha,beta)', ' and diagonals of (S,T)', / '  7 = 1/ULP  if SDIM is not the correct number of ', 'selected eigenvalues', / '  8 = 1/ULP  if DIFEST/DIFTRU > 10*THRESH or ', 'DIFTRU/DIFEST > 10*THRESH', / '  9 = 1/ULP  if DIFEST <> 0 or DIFTRU > ULP*norm(A,B) ', 'when reordering fails', / ' 10 = 1/ULP  if PLEST/PLTRU > THRESH or ', 'PLTRU/PLEST > THRESH', / '    ( Test 10 is only for input examples )', / );
 9991 FORMAT( ' Matrix order=', I2, ', type=', I2, ', a=', D10.3, ', order(A_11)=', I2, ', result ', I2, ' is ', 0P, F8.2 );
 9990 FORMAT( ' Matrix order=', I2, ', type=', I2, ', a=', D10.3, ', order(A_11)=', I2, ', result ', I2, ' is ', 0P, D10.3 );
 9989 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',', ' result ', I2, ' is', 0P, F8.2 );
 9988 FORMAT( ' Input example #', I2, ', matrix order=', I4, ',', ' result ', I2, ' is', 1P, D10.3 );
      }
