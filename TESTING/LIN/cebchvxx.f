      SUBROUTINE CEBCHVXX( THRESH, PATH )
      IMPLICIT NONE
      // .. Scalar Arguments ..
      REAL               THRESH
      String             PATH;

      int                NMAX, NPARAMS, NERRBND, NTESTS, KL, KU;
      const              NMAX = 6, NPARAMS = 2, NERRBND = 3, NTESTS = 6;

      // .. Local Scalars ..
      int                N, NRHS, INFO, I ,J, k, NFAIL, LDA, N_AUX_TESTS, LDAB, LDAFB;
      String             FACT, TRANS, UPLO, EQUED;
      String             C2;
      String   (3)       NGUAR, CGUAR;
      bool               printed_guide;
      REAL               NCOND, CCOND, M, NORMDIF, NORMT, RCOND, RNORM, RINORM, SUMR, SUMRI, EPS, BERR(NMAX), RPVGRW, ORCOND, CWISE_ERR, NWISE_ERR, CWISE_BND, NWISE_BND, CWISE_RCOND, NWISE_RCOND, CONDTHRESH, ERRTHRESH
      COMPLEX            ZDUM

      // .. Local Arrays ..
      REAL               TSTRAT(NTESTS), RINV(NMAX), PARAMS(NPARAMS), S(NMAX), R(NMAX),C(NMAX),RWORK(3*NMAX), DIFF(NMAX, NMAX), ERRBND_N(NMAX*3), ERRBND_C(NMAX*3)
      int                IPIV(NMAX);
      COMPLEX            A(NMAX,NMAX),INVHILB(NMAX,NMAX),X(NMAX,NMAX), WORK(NMAX*3*5), AF(NMAX, NMAX),B(NMAX, NMAX), ACOPY(NMAX, NMAX), AB( (NMAX-1)+(NMAX-1)+1, NMAX ), ABCOPY( (NMAX-1)+(NMAX-1)+1, NMAX ), AFB( 2*(NMAX-1)+(NMAX-1)+1, NMAX )

      // .. External Functions ..
      REAL               SLAMCH

      // .. External Subroutines ..
      // EXTERNAL CLAHILB, CGESVXX, CSYSVXX, CPOSVXX, CGBSVXX, CLACPY, LSAMEN
      bool               LSAMEN;

      // .. Intrinsic Functions ..
      // INTRINSIC SQRT, MAX, ABS, REAL, AIMAG

      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )

      // .. Parameters ..
      int                NWISE_I, CWISE_I;
      const              NWISE_I = 1, CWISE_I = 1;
      int                BND_I, COND_I;
      const              BND_I = 2, COND_I = 3;

*  Create the loop to test out the Hilbert matrices

      FACT = 'E'
      UPLO = 'U'
      TRANS = 'N'
      EQUED = 'N'
      EPS = SLAMCH('Epsilon')
      NFAIL = 0
      N_AUX_TESTS = 0
      LDA = NMAX
      LDAB = (NMAX-1)+(NMAX-1)+1
      LDAFB = 2*(NMAX-1)+(NMAX-1)+1
      C2 = PATH( 2: 3 )

      // Main loop to test the different Hilbert Matrices.

      printed_guide = false;

      for (N = 1; N <= NMAX; N++) {
         PARAMS(1) = -1
         PARAMS(2) = -1

         KL = N-1
         KU = N-1
         NRHS = n
         M = MAX(SQRT(REAL(N)), 10.0)

         // Generate the Hilbert matrix, its inverse, and the
         // right hand side, all scaled by the LCM(1,..,2N-1).
         clahilb(N, N, A, LDA, INVHILB, LDA, B, LDA, WORK, INFO, PATH);

         // Copy A into ACOPY.
         clacpy('ALL', N, N, A, NMAX, ACOPY, NMAX);

         // Store A in band format for GB tests
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= KL+KU+1; I++) {
               AB( I, J ) = (0.0E+0,0.0E+0)
            }
         }
         for (J = 1; J <= N; J++) {
            DO I = MAX( 1, J-KU ), MIN( N, J+KL )
               AB( KU+1+I-J, J ) = A( I, J )
            }
         }

         // Copy AB into ABCOPY.
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= KL+KU+1; I++) {
               ABCOPY( I, J ) = (0.0E+0,0.0E+0)
            }
         }
         clacpy('ALL', KL+KU+1, N, AB, LDAB, ABCOPY, LDAB);

         // Call C**SVXX with default PARAMS and N_ERR_BND = 3.
         if ( LSAMEN( 2, C2, 'SY' ) ) {
            csysvxx(FACT, UPLO, N, NRHS, ACOPY, LDA, AF, LDA, IPIV, EQUED, S, B, LDA, X, LDA, ORCOND, RPVGRW, BERR, NERRBND, ERRBND_N, ERRBND_C, NPARAMS, PARAMS, WORK, RWORK, INFO);
         } else if ( LSAMEN( 2, C2, 'PO' ) ) {
            cposvxx(FACT, UPLO, N, NRHS, ACOPY, LDA, AF, LDA, EQUED, S, B, LDA, X, LDA, ORCOND, RPVGRW, BERR, NERRBND, ERRBND_N, ERRBND_C, NPARAMS, PARAMS, WORK, RWORK, INFO);
         } else if ( LSAMEN( 2, C2, 'HE' ) ) {
            chesvxx(FACT, UPLO, N, NRHS, ACOPY, LDA, AF, LDA, IPIV, EQUED, S, B, LDA, X, LDA, ORCOND, RPVGRW, BERR, NERRBND, ERRBND_N, ERRBND_C, NPARAMS, PARAMS, WORK, RWORK, INFO);
         } else if ( LSAMEN( 2, C2, 'GB' ) ) {
            cgbsvxx(FACT, TRANS, N, KL, KU, NRHS, ABCOPY, LDAB, AFB, LDAFB, IPIV, EQUED, R, C, B, LDA, X, LDA, ORCOND, RPVGRW, BERR, NERRBND, ERRBND_N, ERRBND_C, NPARAMS, PARAMS, WORK, RWORK, INFO);
         } else {
            cgesvxx(FACT, TRANS, N, NRHS, ACOPY, LDA, AF, LDA, IPIV, EQUED, R, C, B, LDA, X, LDA, ORCOND, RPVGRW, BERR, NERRBND, ERRBND_N, ERRBND_C, NPARAMS, PARAMS, WORK, RWORK, INFO);
         }

         N_AUX_TESTS = N_AUX_TESTS + 1
         if (ORCOND .LT. EPS) {
         // Either factorization failed or the matrix is flagged, and 1 <=
         // INFO <= N+1. We don't decide based on rcond anymore.
             // IF (INFO == 0 .OR. INFO .GT. N+1) THEN
                // NFAIL = NFAIL + 1
                // WRITE (*, FMT=8000) N, INFO, ORCOND, RCOND
             // END IF
         } else {
         // Either everything succeeded (INFO == 0) or some solution failed
         // to converge (INFO > N+1).
            if (INFO .GT. 0 .AND. INFO .LE. N+1) {
               NFAIL = NFAIL + 1
               WRITE (*, FMT=8000) C2, N, INFO, ORCOND, RCOND
            }
         }

         // Calculating the difference between C**SVXX's X and the true X.
         for (I = 1; I <= N; I++) {
            for (J = 1; J <= NRHS; J++) {
               DIFF(I,J) = X(I,J) - INVHILB(I,J)
            }
         }

         // Calculating the RCOND
         RNORM = 0
         RINORM = 0
         if ( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'SY' ) .OR. LSAMEN( 2, C2, 'HE' ) ) {
            for (I = 1; I <= N; I++) {
               SUMR = 0
               SUMRI = 0
               for (J = 1; J <= N; J++) {
                  SUMR = SUMR + S(I) * CABS1(A(I,J)) * S(J)
                  SUMRI = SUMRI + CABS1(INVHILB(I, J)) / (S(J) * S(I))
               }
               RNORM = MAX(RNORM,SUMR)
               RINORM = MAX(RINORM,SUMRI)
            }
         } else if ( LSAMEN( 2, C2, 'GE' ) .OR. LSAMEN( 2, C2, 'GB' ) ) {
            for (I = 1; I <= N; I++) {
               SUMR = 0
               SUMRI = 0
               for (J = 1; J <= N; J++) {
                  SUMR = SUMR + R(I) * CABS1(A(I,J)) * C(J)
                  SUMRI = SUMRI + CABS1(INVHILB(I, J)) / (R(J) * C(I))
               }
               RNORM = MAX(RNORM,SUMR)
               RINORM = MAX(RINORM,SUMRI)
            }
         }

         RNORM = RNORM / CABS1(A(1, 1))
         RCOND = 1.0/(RNORM * RINORM)

         // Calculating the R for normwise rcond.
         for (I = 1; I <= N; I++) {
            RINV(I) = 0.0
         }
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= N; I++) {
               RINV(I) = RINV(I) + CABS1(A(I,J))
            }
         }

         // Calculating the Normwise rcond.
         RINORM = 0.0
         for (I = 1; I <= N; I++) {
            SUMRI = 0.0
            for (J = 1; J <= N; J++) {
               SUMRI = SUMRI + CABS1(INVHILB(I,J) * RINV(J))
            }
            RINORM = MAX(RINORM, SUMRI)
         }

         // invhilb is the inverse *unscaled* Hilbert matrix, so scale its norm
         // by 1/A(1,1) to make the scaling match A (the scaled Hilbert matrix)
         NCOND = CABS1(A(1,1)) / RINORM

         CONDTHRESH = M * EPS
         ERRTHRESH = M * EPS

         for (K = 1; K <= NRHS; K++) {
            NORMT = 0.0
            NORMDIF = 0.0
            CWISE_ERR = 0.0
            for (I = 1; I <= N; I++) {
               NORMT = MAX(CABS1(INVHILB(I, K)), NORMT)
               NORMDIF = MAX(CABS1(X(I,K) - INVHILB(I,K)), NORMDIF)
               if (INVHILB(I,K) .NE. 0.0) {
                  CWISE_ERR = MAX(CABS1(X(I,K) - INVHILB(I,K)) /CABS1(INVHILB(I,K)), CWISE_ERR)
               } else if (X(I, K) .NE. 0.0) {
                  CWISE_ERR = SLAMCH('OVERFLOW')
               }
            }
            if (NORMT .NE. 0.0) {
               NWISE_ERR = NORMDIF / NORMT
            } else if (NORMDIF .NE. 0.0) {
               NWISE_ERR = SLAMCH('OVERFLOW')
            } else {
               NWISE_ERR = 0.0
            }

            for (I = 1; I <= N; I++) {
               RINV(I) = 0.0
            }
            for (J = 1; J <= N; J++) {
               for (I = 1; I <= N; I++) {
                  RINV(I) = RINV(I) + CABS1(A(I, J) * INVHILB(J, K))
               }
            }
            RINORM = 0.0
            for (I = 1; I <= N; I++) {
               SUMRI = 0.0
               for (J = 1; J <= N; J++) {
                  SUMRI = SUMRI + CABS1(INVHILB(I, J) * RINV(J) / INVHILB(I, K))
               }
               RINORM = MAX(RINORM, SUMRI)
            }
         // invhilb is the inverse *unscaled* Hilbert matrix, so scale its norm
         // by 1/A(1,1) to make the scaling match A (the scaled Hilbert matrix)
            CCOND = CABS1(A(1,1))/RINORM

         // Forward error bound tests
            NWISE_BND = ERRBND_N(K + (BND_I-1)*NRHS)
            CWISE_BND = ERRBND_C(K + (BND_I-1)*NRHS)
            NWISE_RCOND = ERRBND_N(K + (COND_I-1)*NRHS)
            CWISE_RCOND = ERRBND_C(K + (COND_I-1)*NRHS)
             // write (*,*) 'nwise : ', n, k, ncond, nwise_rcond,
      // $           condthresh, ncond.ge.condthresh
             // write (*,*) 'nwise2: ', k, nwise_bnd, nwise_err, errthresh
            if (NCOND .GE. CONDTHRESH) {
               NGUAR = 'YES'
               if (NWISE_BND .GT. ERRTHRESH) {
                  TSTRAT(1) = 1/(2.0*EPS)
               } else {
                  if (NWISE_BND .NE. 0.0) {
                     TSTRAT(1) = NWISE_ERR / NWISE_BND
                  } else if (NWISE_ERR .NE. 0.0) {
                     TSTRAT(1) = 1/(16.0*EPS)
                  } else {
                     TSTRAT(1) = 0.0
                  }
                  if (TSTRAT(1) .GT. 1.0) {
                     TSTRAT(1) = 1/(4.0*EPS)
                  }
               }
            } else {
               NGUAR = 'NO'
               if (NWISE_BND .LT. 1.0) {
                  TSTRAT(1) = 1/(8.0*EPS)
               } else {
                  TSTRAT(1) = 1.0
               }
            }
             // write (*,*) 'cwise : ', n, k, ccond, cwise_rcond,
      // $           condthresh, ccond.ge.condthresh
             // write (*,*) 'cwise2: ', k, cwise_bnd, cwise_err, errthresh
            if (CCOND .GE. CONDTHRESH) {
               CGUAR = 'YES'
               if (CWISE_BND .GT. ERRTHRESH) {
                  TSTRAT(2) = 1/(2.0*EPS)
               } else {
                  if (CWISE_BND .NE. 0.0) {
                     TSTRAT(2) = CWISE_ERR / CWISE_BND
                  } else if (CWISE_ERR .NE. 0.0) {
                     TSTRAT(2) = 1/(16.0*EPS)
                  } else {
                     TSTRAT(2) = 0.0
                  }
                  IF (TSTRAT(2) .GT. 1.0) TSTRAT(2) = 1/(4.0*EPS)
               }
            } else {
               CGUAR = 'NO'
               if (CWISE_BND .LT. 1.0) {
                  TSTRAT(2) = 1/(8.0*EPS)
               } else {
                  TSTRAT(2) = 1.0
               }
            }

      // Backwards error test
            TSTRAT(3) = BERR(K)/EPS

      // Condition number tests
            TSTRAT(4) = RCOND / ORCOND
            IF (RCOND .GE. CONDTHRESH .AND. TSTRAT(4) .LT. 1.0) TSTRAT(4) = 1.0 / TSTRAT(4)

            TSTRAT(5) = NCOND / NWISE_RCOND
            IF (NCOND .GE. CONDTHRESH .AND. TSTRAT(5) .LT. 1.0) TSTRAT(5) = 1.0 / TSTRAT(5)

            TSTRAT(6) = CCOND / NWISE_RCOND
            IF (CCOND .GE. CONDTHRESH .AND. TSTRAT(6) .LT. 1.0) TSTRAT(6) = 1.0 / TSTRAT(6)

            for (I = 1; I <= NTESTS; I++) {
               if (TSTRAT(I) .GT. THRESH) {
                  if (.NOT.PRINTED_GUIDE) {
                     WRITE(*,*)
                     WRITE( *, 9996) 1
                     WRITE( *, 9995) 2
                     WRITE( *, 9994) 3
                     WRITE( *, 9993) 4
                     WRITE( *, 9992) 5
                     WRITE( *, 9991) 6
                     WRITE( *, 9990) 7
                     WRITE( *, 9989) 8
                     WRITE(*,*)
                     PRINTED_GUIDE = true;
                  }
                  WRITE( *, 9999) C2, N, K, NGUAR, CGUAR, I, TSTRAT(I)
                  NFAIL = NFAIL + 1
               }
            }
      }

c$$$         WRITE(*,*)
c$$$         WRITE(*,*) 'Normwise Error Bounds'
c$$$         WRITE(*,*) 'Guaranteed error bound: ',ERRBND(NRHS,nwise_i,bnd_i)
c$$$         WRITE(*,*) 'Reciprocal condition number: ',ERRBND(NRHS,nwise_i,cond_i)
c$$$         WRITE(*,*) 'Raw error estimate: ',ERRBND(NRHS,nwise_i,rawbnd_i)
c$$$         WRITE(*,*)
c$$$         WRITE(*,*) 'Componentwise Error Bounds'
c$$$         WRITE(*,*) 'Guaranteed error bound: ',ERRBND(NRHS,cwise_i,bnd_i)
c$$$         WRITE(*,*) 'Reciprocal condition number: ',ERRBND(NRHS,cwise_i,cond_i)
c$$$         WRITE(*,*) 'Raw error estimate: ',ERRBND(NRHS,cwise_i,rawbnd_i)
c$$$         print *, 'Info: ', info
c$$$         WRITE(*,*)
          // WRITE(*,*) 'TSTRAT: ',TSTRAT

      }

      WRITE(*,*)
      if ( NFAIL .GT. 0 ) {
         WRITE(*,9998) C2, NFAIL, NTESTS*N+N_AUX_TESTS
      } else {
         WRITE(*,9997) C2
      }
 9999 FORMAT( ' C', A2, 'SVXX: N =', I2, ', RHS = ', I2, ', NWISE GUAR. = ', A, ', CWISE GUAR. = ', A, ' test(',I1,') =', G12.5 )
 9998 FORMAT( ' C', A2, 'SVXX: ', I6, ' out of ', I6, ' tests failed to pass the threshold' )
 9997 FORMAT( ' C', A2, 'SVXX passed the tests of error bounds' )
      // Test ratios.
 9996 FORMAT( 3X, I2, ': Normwise guaranteed forward error', / 5X, 'Guaranteed case: if norm ( abs( Xc - Xt )', ' / norm ( Xt ) .LE. ERRBND( *, nwise_i, bnd_i ), then', / 5X, 'ERRBND( *, nwise_i, bnd_i ) .LE. MAX(SQRT(N), 10) * EPS')
 9995 FORMAT( 3X, I2, ': Componentwise guaranteed forward error' )
 9994 FORMAT( 3X, I2, ': Backwards error' )
 9993 FORMAT( 3X, I2, ': Reciprocal condition number' )
 9992 FORMAT( 3X, I2, ': Reciprocal normwise condition number' )
 9991 FORMAT( 3X, I2, ': Raw normwise error estimate' )
 9990 FORMAT( 3X, I2, ': Reciprocal componentwise condition number' )
 9989 FORMAT( 3X, I2, ': Raw componentwise error estimate' )

 8000 FORMAT( ' C', A2, 'SVXX: N =', I2, ', INFO = ', I3, ', ORCOND = ', G12.5, ', real RCOND = ', G12.5 )

      // End of CEBCHVXX

      }
