      SUBROUTINE CLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, RANK, WORK, RWORK, IWORK, INFO )

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
      REAL               D( * ), E( * ), RWORK( * )
      COMPLEX            B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 ;
      COMPLEX            CZERO
      const              CZERO = ( 0.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      int                BX, BXST, C, DIFL, DIFR, GIVCOL, GIVNUM, GIVPTR, I, ICMPQ1, ICMPQ2, IRWB, IRWIB, IRWRB, IRWU, IRWVT, IRWWRK, IWK, J, JCOL, JIMAG, JREAL, JROW, K, NLVL, NM1, NRWORK, NSIZE, NSUB, PERM, POLES, S, SIZEI, SMLSZP, SQRE, ST, ST1, U, VT, Z;
      REAL               CS, EPS, ORGNRM, R, RCND, SN, TOL
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SLAMCH, SLANST
      // EXTERNAL ISAMAX, SLAMCH, SLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACPY, CLALSA, CLASCL, CLASET, CSROT, SGEMM, SLARTG, SLASCL, SLASDA, SLASDQ, SLASET, SLASRT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, INT, LOG, REAL, SIGN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( N.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.1 ) {
         INFO = -4
      } else if ( ( LDB.LT.1 ) .OR. ( LDB.LT.N ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('CLALSD', -INFO );
         RETURN
      }

      EPS = SLAMCH( 'Epsilon' )

      // Set up the tolerance.

      if ( ( RCOND.LE.ZERO ) .OR. ( RCOND.GE.ONE ) ) {
         RCND = EPS
      } else {
         RCND = RCOND
      }

      RANK = 0

      // Quick return if possible.

      if ( N.EQ.0 ) {
         RETURN
      } else if ( N.EQ.1 ) {
         if ( D( 1 ).EQ.ZERO ) {
            claset('A', 1, NRHS, CZERO, CZERO, B, LDB );
         } else {
            RANK = 1
            clascl('G', 0, 0, D( 1 ), ONE, 1, NRHS, B, LDB, INFO );
            D( 1 ) = ABS( D( 1 ) )
         }
         RETURN
      }

      // Rotate the matrix if it is lower bidiagonal.

      if ( UPLO.EQ.'L' ) {
         for (I = 1; I <= N - 1; I++) { // 10
            slartg(D( I ), E( I ), CS, SN, R );
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            if ( NRHS.EQ.1 ) {
               csrot(1, B( I, 1 ), 1, B( I+1, 1 ), 1, CS, SN );
            } else {
               RWORK( I*2-1 ) = CS
               RWORK( I*2 ) = SN
            }
         } // 10
         if ( NRHS.GT.1 ) {
            for (I = 1; I <= NRHS; I++) { // 30
               for (J = 1; J <= N - 1; J++) { // 20
                  CS = RWORK( J*2-1 )
                  SN = RWORK( J*2 )
                  csrot(1, B( J, I ), 1, B( J+1, I ), 1, CS, SN );
               } // 20
            } // 30
         }
      }

      // Scale.

      NM1 = N - 1
      ORGNRM = SLANST( 'M', N, D, E )
      if ( ORGNRM.EQ.ZERO ) {
         claset('A', N, NRHS, CZERO, CZERO, B, LDB );
         RETURN
      }

      slascl('G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO );
      slascl('G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, INFO );

      // If N is smaller than the minimum divide size SMLSIZ, then solve
      // the problem with another solver.

      if ( N.LE.SMLSIZ ) {
         IRWU = 1
         IRWVT = IRWU + N*N
         IRWWRK = IRWVT + N*N
         IRWRB = IRWWRK
         IRWIB = IRWRB + N*NRHS
         IRWB = IRWIB + N*NRHS
         slaset('A', N, N, ZERO, ONE, RWORK( IRWU ), N );
         slaset('A', N, N, ZERO, ONE, RWORK( IRWVT ), N );
         slasdq('U', 0, N, N, N, 0, D, E, RWORK( IRWVT ), N, RWORK( IRWU ), N, RWORK( IRWWRK ), 1, RWORK( IRWWRK ), INFO );
         if ( INFO.NE.0 ) {
            RETURN
         }

         // In the real version, B is passed to SLASDQ and multiplied
         // internally by Q**H. Here B is complex and that product is
         // computed below in two steps (real and imaginary parts).

         J = IRWB - 1
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 50
            for (JROW = 1; JROW <= N; JROW++) { // 40
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
            } // 40
         } // 50
         sgemm('T', 'N', N, NRHS, N, ONE, RWORK( IRWU ), N, RWORK( IRWB ), N, ZERO, RWORK( IRWRB ), N );
         J = IRWB - 1
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 70
            for (JROW = 1; JROW <= N; JROW++) { // 60
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
            } // 60
         } // 70
         sgemm('T', 'N', N, NRHS, N, ONE, RWORK( IRWU ), N, RWORK( IRWB ), N, ZERO, RWORK( IRWIB ), N );
         JREAL = IRWRB - 1
         JIMAG = IRWIB - 1
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 90
            for (JROW = 1; JROW <= N; JROW++) { // 80
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               B( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
            } // 80
         } // 90

         TOL = RCND*ABS( D( ISAMAX( N, D, 1 ) ) )
         for (I = 1; I <= N; I++) { // 100
            if ( D( I ).LE.TOL ) {
               claset('A', 1, NRHS, CZERO, CZERO, B( I, 1 ), LDB );
            } else {
               clascl('G', 0, 0, D( I ), ONE, 1, NRHS, B( I, 1 ), LDB, INFO );
               RANK = RANK + 1
            }
         } // 100

         // Since B is complex, the following call to SGEMM is performed
         // in two steps (real and imaginary parts). That is for V * B
         // (in the real version of the code V**H is stored in WORK).

         // CALL SGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO,
*    $               WORK( NWORK ), N )

         J = IRWB - 1
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 120
            for (JROW = 1; JROW <= N; JROW++) { // 110
               J = J + 1
               RWORK( J ) = REAL( B( JROW, JCOL ) )
            } // 110
         } // 120
         sgemm('T', 'N', N, NRHS, N, ONE, RWORK( IRWVT ), N, RWORK( IRWB ), N, ZERO, RWORK( IRWRB ), N );
         J = IRWB - 1
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 140
            for (JROW = 1; JROW <= N; JROW++) { // 130
               J = J + 1
               RWORK( J ) = AIMAG( B( JROW, JCOL ) )
            } // 130
         } // 140
         sgemm('T', 'N', N, NRHS, N, ONE, RWORK( IRWVT ), N, RWORK( IRWB ), N, ZERO, RWORK( IRWIB ), N );
         JREAL = IRWRB - 1
         JIMAG = IRWIB - 1
         for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 160
            for (JROW = 1; JROW <= N; JROW++) { // 150
               JREAL = JREAL + 1
               JIMAG = JIMAG + 1
               B( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
            } // 150
         } // 160

         // Unscale.

         slascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO );
         slasrt('D', N, D, INFO );
         clascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO );

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
      NRWORK = GIVNUM + 2*NLVL*N
      BX = 1

      IRWRB = NRWORK
      IRWIB = IRWRB + SMLSIZ*NRHS
      IRWB = IRWIB + SMLSIZ*NRHS

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

      for (I = 1; I <= N; I++) { // 170
         if ( ABS( D( I ) ).LT.EPS ) {
            D( I ) = SIGN( EPS, D( I ) )
         }
      } // 170

      for (I = 1; I <= NM1; I++) { // 240
         if ( ( ABS( E( I ) ).LT.EPS ) .OR. ( I.EQ.NM1 ) ) {
            NSUB = NSUB + 1
            IWORK( NSUB ) = ST

            // Subproblem found. First determine its size and then
            // apply divide and conquer on it.

            if ( I.LT.NM1 ) {

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
               ccopy(NRHS, B( N, 1 ), LDB, WORK( BX+NM1 ), N );
            }
            ST1 = ST - 1
            if ( NSIZE.EQ.1 ) {

               // This is a 1-by-1 subproblem and is not solved
               // explicitly.

               ccopy(NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N );
            } else if ( NSIZE.LE.SMLSIZ ) {

               // This is a small subproblem and is solved by SLASDQ.

               slaset('A', NSIZE, NSIZE, ZERO, ONE, RWORK( VT+ST1 ), N );
               slaset('A', NSIZE, NSIZE, ZERO, ONE, RWORK( U+ST1 ), N );
               slasdq('U', 0, NSIZE, NSIZE, NSIZE, 0, D( ST ), E( ST ), RWORK( VT+ST1 ), N, RWORK( U+ST1 ), N, RWORK( NRWORK ), 1, RWORK( NRWORK ), INFO );
               if ( INFO.NE.0 ) {
                  RETURN
               }

               // In the real version, B is passed to SLASDQ and multiplied
               // internally by Q**H. Here B is complex and that product is
               // computed below in two steps (real and imaginary parts).

               J = IRWB - 1
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 190
                  for (JROW = ST; JROW <= ST + NSIZE - 1; JROW++) { // 180
                     J = J + 1
                     RWORK( J ) = REAL( B( JROW, JCOL ) )
                  } // 180
               } // 190
               sgemm('T', 'N', NSIZE, NRHS, NSIZE, ONE, RWORK( U+ST1 ), N, RWORK( IRWB ), NSIZE, ZERO, RWORK( IRWRB ), NSIZE );
               J = IRWB - 1
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 210
                  for (JROW = ST; JROW <= ST + NSIZE - 1; JROW++) { // 200
                     J = J + 1
                     RWORK( J ) = AIMAG( B( JROW, JCOL ) )
                  } // 200
               } // 210
               sgemm('T', 'N', NSIZE, NRHS, NSIZE, ONE, RWORK( U+ST1 ), N, RWORK( IRWB ), NSIZE, ZERO, RWORK( IRWIB ), NSIZE );
               JREAL = IRWRB - 1
               JIMAG = IRWIB - 1
               for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 230
                  for (JROW = ST; JROW <= ST + NSIZE - 1; JROW++) { // 220
                     JREAL = JREAL + 1
                     JIMAG = JIMAG + 1
                     B( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
                  } // 220
               } // 230

               clacpy('A', NSIZE, NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N );
            } else {

               // A large problem. Solve it using divide and conquer.

               slasda(ICMPQ1, SMLSIZ, NSIZE, SQRE, D( ST ), E( ST ), RWORK( U+ST1 ), N, RWORK( VT+ST1 ), IWORK( K+ST1 ), RWORK( DIFL+ST1 ), RWORK( DIFR+ST1 ), RWORK( Z+ST1 ), RWORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), RWORK( GIVNUM+ST1 ), RWORK( C+ST1 ), RWORK( S+ST1 ), RWORK( NRWORK ), IWORK( IWK ), INFO );
               if ( INFO.NE.0 ) {
                  RETURN
               }
               BXST = BX + ST1
               clalsa(ICMPQ2, SMLSIZ, NSIZE, NRHS, B( ST, 1 ), LDB, WORK( BXST ), N, RWORK( U+ST1 ), N, RWORK( VT+ST1 ), IWORK( K+ST1 ), RWORK( DIFL+ST1 ), RWORK( DIFR+ST1 ), RWORK( Z+ST1 ), RWORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), RWORK( GIVNUM+ST1 ), RWORK( C+ST1 ), RWORK( S+ST1 ), RWORK( NRWORK ), IWORK( IWK ), INFO );
               if ( INFO.NE.0 ) {
                  RETURN
               }
            }
            ST = I + 1
         }
      } // 240

      // Apply the singular values and treat the tiny ones as zero.

      TOL = RCND*ABS( D( ISAMAX( N, D, 1 ) ) )

      for (I = 1; I <= N; I++) { // 250

         // Some of the elements in D can be negative because 1-by-1
         // subproblems were not solved explicitly.

         if ( ABS( D( I ) ).LE.TOL ) {
            claset('A', 1, NRHS, CZERO, CZERO, WORK( BX+I-1 ), N );
         } else {
            RANK = RANK + 1
            clascl('G', 0, 0, D( I ), ONE, 1, NRHS, WORK( BX+I-1 ), N, INFO );
         }
         D( I ) = ABS( D( I ) )
      } // 250

      // Now apply back the right singular vectors.

      ICMPQ2 = 1
      for (I = 1; I <= NSUB; I++) { // 320
         ST = IWORK( I )
         ST1 = ST - 1
         NSIZE = IWORK( SIZEI+I-1 )
         BXST = BX + ST1
         if ( NSIZE.EQ.1 ) {
            ccopy(NRHS, WORK( BXST ), N, B( ST, 1 ), LDB );
         } else if ( NSIZE.LE.SMLSIZ ) {

            // Since B and BX are complex, the following call to SGEMM
            // is performed in two steps (real and imaginary parts).

            // CALL SGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE,
*    $                  RWORK( VT+ST1 ), N, RWORK( BXST ), N, ZERO,
*    $                  B( ST, 1 ), LDB )

            J = BXST - N - 1
            JREAL = IRWB - 1
            for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 270
               J = J + N
               for (JROW = 1; JROW <= NSIZE; JROW++) { // 260
                  JREAL = JREAL + 1
                  RWORK( JREAL ) = REAL( WORK( J+JROW ) )
               } // 260
            } // 270
            sgemm('T', 'N', NSIZE, NRHS, NSIZE, ONE, RWORK( VT+ST1 ), N, RWORK( IRWB ), NSIZE, ZERO, RWORK( IRWRB ), NSIZE );
            J = BXST - N - 1
            JIMAG = IRWB - 1
            for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 290
               J = J + N
               for (JROW = 1; JROW <= NSIZE; JROW++) { // 280
                  JIMAG = JIMAG + 1
                  RWORK( JIMAG ) = AIMAG( WORK( J+JROW ) )
               } // 280
            } // 290
            sgemm('T', 'N', NSIZE, NRHS, NSIZE, ONE, RWORK( VT+ST1 ), N, RWORK( IRWB ), NSIZE, ZERO, RWORK( IRWIB ), NSIZE );
            JREAL = IRWRB - 1
            JIMAG = IRWIB - 1
            for (JCOL = 1; JCOL <= NRHS; JCOL++) { // 310
               for (JROW = ST; JROW <= ST + NSIZE - 1; JROW++) { // 300
                  JREAL = JREAL + 1
                  JIMAG = JIMAG + 1
                  B( JROW, JCOL ) = CMPLX( RWORK( JREAL ), RWORK( JIMAG ) )
               } // 300
            } // 310
         } else {
            clalsa(ICMPQ2, SMLSIZ, NSIZE, NRHS, WORK( BXST ), N, B( ST, 1 ), LDB, RWORK( U+ST1 ), N, RWORK( VT+ST1 ), IWORK( K+ST1 ), RWORK( DIFL+ST1 ), RWORK( DIFR+ST1 ), RWORK( Z+ST1 ), RWORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), RWORK( GIVNUM+ST1 ), RWORK( C+ST1 ), RWORK( S+ST1 ), RWORK( NRWORK ), IWORK( IWK ), INFO );
            if ( INFO.NE.0 ) {
               RETURN
            }
         }
      } // 320

      // Unscale and sort the singular values.

      slascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO );
      slasrt('D', N, D, INFO );
      clascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO );

      RETURN

      // End of CLALSD

      }
