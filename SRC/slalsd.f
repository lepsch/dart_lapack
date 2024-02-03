      SUBROUTINE SLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, RANK, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS, RANK, SMLSIZ;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               B( LDB, * ), D( * ), E( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 ;
      // ..
      // .. Local Scalars ..
      int                BX, BXST, C, DIFL, DIFR, GIVCOL, GIVNUM, GIVPTR, I, ICMPQ1, ICMPQ2, IWK, J, K, NLVL, NM1, NSIZE, NSUB, NWORK, PERM, POLES, S, SIZEI, SMLSZP, SQRE, ST, ST1, U, VT, Z;
      REAL               CS, EPS, ORGNRM, R, RCND, SN, TOL
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SLAMCH, SLANST
      // EXTERNAL ISAMAX, SLAMCH, SLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMM, SLACPY, SLALSA, SLARTG, SLASCL, SLASDA, SLASDQ, SLASET, SLASRT, SROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, REAL, SIGN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( N < 0 ) {
         INFO = -3
      } else if ( NRHS < 1 ) {
         INFO = -4
      } else if ( ( LDB < 1 ) || ( LDB < N ) ) {
         INFO = -8
      }
      if ( INFO != 0 ) {
         xerbla('SLALSD', -INFO );
         RETURN
      }

      EPS = SLAMCH( 'Epsilon' )

      // Set up the tolerance.

      if ( ( RCOND.LE.ZERO ) || ( RCOND.GE.ONE ) ) {
         RCND = EPS
      } else {
         RCND = RCOND
      }

      RANK = 0

      // Quick return if possible.

      if ( N == 0 ) {
         RETURN
      } else if ( N == 1 ) {
         if ( D( 1 ) == ZERO ) {
            slaset('A', 1, NRHS, ZERO, ZERO, B, LDB );
         } else {
            RANK = 1
            slascl('G', 0, 0, D( 1 ), ONE, 1, NRHS, B, LDB, INFO );
            D( 1 ) = ABS( D( 1 ) )
         }
         RETURN
      }

      // Rotate the matrix if it is lower bidiagonal.

      if ( UPLO == 'L' ) {
         for (I = 1; I <= N - 1; I++) { // 10
            slartg(D( I ), E( I ), CS, SN, R );
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            if ( NRHS == 1 ) {
               srot(1, B( I, 1 ), 1, B( I+1, 1 ), 1, CS, SN );
            } else {
               WORK( I*2-1 ) = CS
               WORK( I*2 ) = SN
            }
         } // 10
         if ( NRHS.GT.1 ) {
            for (I = 1; I <= NRHS; I++) { // 30
               for (J = 1; J <= N - 1; J++) { // 20
                  CS = WORK( J*2-1 )
                  SN = WORK( J*2 )
                  srot(1, B( J, I ), 1, B( J+1, I ), 1, CS, SN );
               } // 20
            } // 30
         }
      }

      // Scale.

      NM1 = N - 1
      ORGNRM = SLANST( 'M', N, D, E )
      if ( ORGNRM == ZERO ) {
         slaset('A', N, NRHS, ZERO, ZERO, B, LDB );
         RETURN
      }

      slascl('G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO );
      slascl('G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, INFO );

      // If N is smaller than the minimum divide size SMLSIZ, then solve
      // the problem with another solver.

      if ( N.LE.SMLSIZ ) {
         NWORK = 1 + N*N
         slaset('A', N, N, ZERO, ONE, WORK, N );
         slasdq('U', 0, N, N, 0, NRHS, D, E, WORK, N, WORK, N, B, LDB, WORK( NWORK ), INFO );
         if ( INFO != 0 ) {
            RETURN
         }
         TOL = RCND*ABS( D( ISAMAX( N, D, 1 ) ) )
         for (I = 1; I <= N; I++) { // 40
            if ( D( I ).LE.TOL ) {
               slaset('A', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB );
            } else {
               slascl('G', 0, 0, D( I ), ONE, 1, NRHS, B( I, 1 ), LDB, INFO );
               RANK = RANK + 1
            }
         } // 40
         sgemm('T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO, WORK( NWORK ), N );
         slacpy('A', N, NRHS, WORK( NWORK ), N, B, LDB );

         // Unscale.

         slascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO );
         slasrt('D', N, D, INFO );
         slascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO );

         RETURN
      }

      // Book-keeping and setting up some constants.

      NLVL = INT( LOG( REAL( N ) / REAL( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1

      SMLSZP = SMLSIZ + 1

      U = 1
      VT = 1 + SMLSIZ*N
      DIFL = VT + SMLSZP*N
      DIFR = DIFL + NLVL*N
      Z = DIFR + NLVL*N*2
      C = Z + NLVL*N
      S = C + N
      POLES = S + N
      GIVNUM = POLES + 2*NLVL*N
      BX = GIVNUM + 2*NLVL*N
      NWORK = BX + N*NRHS

      SIZEI = 1 + N
      K = SIZEI + N
      GIVPTR = K + N
      PERM = GIVPTR + N
      GIVCOL = PERM + NLVL*N
      IWK = GIVCOL + NLVL*N*2

      ST = 1
      SQRE = 0
      ICMPQ1 = 1
      ICMPQ2 = 0
      NSUB = 0

      for (I = 1; I <= N; I++) { // 50
         if ( ABS( D( I ) ) < EPS ) {
            D( I ) = SIGN( EPS, D( I ) )
         }
      } // 50

      for (I = 1; I <= NM1; I++) { // 60
         if ( ( ABS( E( I ) ) < EPS ) || ( I == NM1 ) ) {
            NSUB = NSUB + 1
            IWORK( NSUB ) = ST

            // Subproblem found. First determine its size and then
            // apply divide and conquer on it.

            if ( I < NM1 ) {

               // A subproblem with E(I) small for I < NM1.

               NSIZE = I - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
            } else if ( ABS( E( I ) ).GE.EPS ) {

               // A subproblem with E(NM1) not too small but I = NM1.

               NSIZE = N - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
            } else {

               // A subproblem with E(NM1) small. This implies an
               // 1-by-1 subproblem at D(N), which is not solved
               // explicitly.

               NSIZE = I - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
               NSUB = NSUB + 1
               IWORK( NSUB ) = N
               IWORK( SIZEI+NSUB-1 ) = 1
               scopy(NRHS, B( N, 1 ), LDB, WORK( BX+NM1 ), N );
            }
            ST1 = ST - 1
            if ( NSIZE == 1 ) {

               // This is a 1-by-1 subproblem and is not solved
               // explicitly.

               scopy(NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N );
            } else if ( NSIZE.LE.SMLSIZ ) {

               // This is a small subproblem and is solved by SLASDQ.

               slaset('A', NSIZE, NSIZE, ZERO, ONE, WORK( VT+ST1 ), N );
               slasdq('U', 0, NSIZE, NSIZE, 0, NRHS, D( ST ), E( ST ), WORK( VT+ST1 ), N, WORK( NWORK ), N, B( ST, 1 ), LDB, WORK( NWORK ), INFO );
               if ( INFO != 0 ) {
                  RETURN
               }
               slacpy('A', NSIZE, NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N );
            } else {

               // A large problem. Solve it using divide and conquer.

               slasda(ICMPQ1, SMLSIZ, NSIZE, SQRE, D( ST ), E( ST ), WORK( U+ST1 ), N, WORK( VT+ST1 ), IWORK( K+ST1 ), WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), WORK( Z+ST1 ), WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ), INFO );
               if ( INFO != 0 ) {
                  RETURN
               }
               BXST = BX + ST1
               slalsa(ICMPQ2, SMLSIZ, NSIZE, NRHS, B( ST, 1 ), LDB, WORK( BXST ), N, WORK( U+ST1 ), N, WORK( VT+ST1 ), IWORK( K+ST1 ), WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), WORK( Z+ST1 ), WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ), INFO );
               if ( INFO != 0 ) {
                  RETURN
               }
            }
            ST = I + 1
         }
      } // 60

      // Apply the singular values and treat the tiny ones as zero.

      TOL = RCND*ABS( D( ISAMAX( N, D, 1 ) ) )

      for (I = 1; I <= N; I++) { // 70

         // Some of the elements in D can be negative because 1-by-1
         // subproblems were not solved explicitly.

         if ( ABS( D( I ) ).LE.TOL ) {
            slaset('A', 1, NRHS, ZERO, ZERO, WORK( BX+I-1 ), N );
         } else {
            RANK = RANK + 1
            slascl('G', 0, 0, D( I ), ONE, 1, NRHS, WORK( BX+I-1 ), N, INFO );
         }
         D( I ) = ABS( D( I ) )
      } // 70

      // Now apply back the right singular vectors.

      ICMPQ2 = 1
      for (I = 1; I <= NSUB; I++) { // 80
         ST = IWORK( I )
         ST1 = ST - 1
         NSIZE = IWORK( SIZEI+I-1 )
         BXST = BX + ST1
         if ( NSIZE == 1 ) {
            scopy(NRHS, WORK( BXST ), N, B( ST, 1 ), LDB );
         } else if ( NSIZE.LE.SMLSIZ ) {
            sgemm('T', 'N', NSIZE, NRHS, NSIZE, ONE, WORK( VT+ST1 ), N, WORK( BXST ), N, ZERO, B( ST, 1 ), LDB );
         } else {
            slalsa(ICMPQ2, SMLSIZ, NSIZE, NRHS, WORK( BXST ), N, B( ST, 1 ), LDB, WORK( U+ST1 ), N, WORK( VT+ST1 ), IWORK( K+ST1 ), WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), WORK( Z+ST1 ), WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ), INFO );
            if ( INFO != 0 ) {
               RETURN
            }
         }
      } // 80

      // Unscale and sort the singular values.

      slascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO );
      slasrt('D', N, D, INFO );
      slascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO );

      RETURN

      // End of SLALSD

      }
