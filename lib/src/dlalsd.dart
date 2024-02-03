      void dlalsd(UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, RANK, WORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS, RANK, SMLSIZ;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             B( LDB, * ), D( * ), E( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      int                BX, BXST, C, DIFL, DIFR, GIVCOL, GIVNUM, GIVPTR, I, ICMPQ1, ICMPQ2, IWK, J, K, NLVL, NM1, NSIZE, NSUB, NWORK, PERM, POLES, S, SIZEI, SMLSZP, SQRE, ST, ST1, U, VT, Z;
      double             CS, EPS, ORGNRM, R, RCND, SN, TOL;
      // ..
      // .. External Functions ..
      //- int                IDAMAX;
      //- double             DLAMCH, DLANST;
      // EXTERNAL IDAMAX, DLAMCH, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DLACPY, DLALSA, DLARTG, DLASCL, DLASDA, DLASDQ, DLASET, DLASRT, DROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, LOG, SIGN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 1 ) {
         INFO = -4;
      } else if ( ( LDB < 1 ) || ( LDB < N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('DLALSD', -INFO );
         return;
      }

      EPS = DLAMCH( 'Epsilon' );

      // Set up the tolerance.

      if ( ( RCOND <= ZERO ) || ( RCOND >= ONE ) ) {
         RCND = EPS;
      } else {
         RCND = RCOND;
      }

      RANK = 0;

      // Quick return if possible.

      if ( N == 0 ) {
         return;
      } else if ( N == 1 ) {
         if ( D( 1 ) == ZERO ) {
            dlaset('A', 1, NRHS, ZERO, ZERO, B, LDB );
         } else {
            RANK = 1;
            dlascl('G', 0, 0, D( 1 ), ONE, 1, NRHS, B, LDB, INFO );
            D( 1 ) = ( D( 1 ) ).abs();
         }
         return;
      }

      // Rotate the matrix if it is lower bidiagonal.

      if ( UPLO == 'L' ) {
         for (I = 1; I <= N - 1; I++) { // 10
            dlartg(D( I ), E( I ), CS, SN, R );
            D( I ) = R;
            E( I ) = SN*D( I+1 );
            D( I+1 ) = CS*D( I+1 );
            if ( NRHS == 1 ) {
               drot(1, B( I, 1 ), 1, B( I+1, 1 ), 1, CS, SN );
            } else {
               WORK( I*2-1 ) = CS;
               WORK( I*2 ) = SN;
            }
         } // 10
         if ( NRHS > 1 ) {
            for (I = 1; I <= NRHS; I++) { // 30
               for (J = 1; J <= N - 1; J++) { // 20
                  CS = WORK( J*2-1 );
                  SN = WORK( J*2 );
                  drot(1, B( J, I ), 1, B( J+1, I ), 1, CS, SN );
               } // 20
            } // 30
         }
      }

      // Scale.

      NM1 = N - 1;
      ORGNRM = DLANST( 'M', N, D, E );
      if ( ORGNRM == ZERO ) {
         dlaset('A', N, NRHS, ZERO, ZERO, B, LDB );
         return;
      }

      dlascl('G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO );
      dlascl('G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, INFO );

      // If N is smaller than the minimum divide size SMLSIZ, then solve
      // the problem with another solver.

      if ( N <= SMLSIZ ) {
         NWORK = 1 + N*N;
         dlaset('A', N, N, ZERO, ONE, WORK, N );
         dlasdq('U', 0, N, N, 0, NRHS, D, E, WORK, N, WORK, N, B, LDB, WORK( NWORK ), INFO );
         if ( INFO != 0 ) {
            return;
         }
         TOL = RCND*ABS( D( IDAMAX( N, D, 1 ) ) );
         for (I = 1; I <= N; I++) { // 40
            if ( D( I ) <= TOL ) {
               dlaset('A', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB );
            } else {
               dlascl('G', 0, 0, D( I ), ONE, 1, NRHS, B( I, 1 ), LDB, INFO );
               RANK = RANK + 1;
            }
         } // 40
         dgemm('T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO, WORK( NWORK ), N );
         dlacpy('A', N, NRHS, WORK( NWORK ), N, B, LDB );

         // Unscale.

         dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO );
         dlasrt('D', N, D, INFO );
         dlascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO );

         return;
      }

      // Book-keeping and setting up some constants.

      NLVL = INT( LOG( DBLE( N ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1;

      SMLSZP = SMLSIZ + 1;

      U = 1;
      VT = 1 + SMLSIZ*N;
      DIFL = VT + SMLSZP*N;
      DIFR = DIFL + NLVL*N;
      Z = DIFR + NLVL*N*2;
      C = Z + NLVL*N;
      S = C + N;
      POLES = S + N;
      GIVNUM = POLES + 2*NLVL*N;
      BX = GIVNUM + 2*NLVL*N;
      NWORK = BX + N*NRHS;

      SIZEI = 1 + N;
      K = SIZEI + N;
      GIVPTR = K + N;
      PERM = GIVPTR + N;
      GIVCOL = PERM + NLVL*N;
      IWK = GIVCOL + NLVL*N*2;

      ST = 1;
      SQRE = 0;
      ICMPQ1 = 1;
      ICMPQ2 = 0;
      NSUB = 0;

      for (I = 1; I <= N; I++) { // 50
         if ( ( D( I ) ).abs() < EPS ) {
            D( I ) = SIGN( EPS, D( I ) );
         }
      } // 50

      for (I = 1; I <= NM1; I++) { // 60
         if ( ( ( E( I ) ).abs() < EPS ) || ( I == NM1 ) ) {
            NSUB = NSUB + 1;
            IWORK( NSUB ) = ST;

            // Subproblem found. First determine its size and then
            // apply divide and conquer on it.

            if ( I < NM1 ) {

               // A subproblem with E(I) small for I < NM1.

               NSIZE = I - ST + 1;
               IWORK( SIZEI+NSUB-1 ) = NSIZE;
            } else if ( ( E( I ) ).abs() >= EPS ) {

               // A subproblem with E(NM1) not too small but I = NM1.

               NSIZE = N - ST + 1;
               IWORK( SIZEI+NSUB-1 ) = NSIZE;
            } else {

               // A subproblem with E(NM1) small. This implies an
               // 1-by-1 subproblem at D(N), which is not solved
               // explicitly.

               NSIZE = I - ST + 1;
               IWORK( SIZEI+NSUB-1 ) = NSIZE;
               NSUB = NSUB + 1;
               IWORK( NSUB ) = N;
               IWORK( SIZEI+NSUB-1 ) = 1;
               dcopy(NRHS, B( N, 1 ), LDB, WORK( BX+NM1 ), N );
            }
            ST1 = ST - 1;
            if ( NSIZE == 1 ) {

               // This is a 1-by-1 subproblem and is not solved
               // explicitly.

               dcopy(NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N );
            } else if ( NSIZE <= SMLSIZ ) {

               // This is a small subproblem and is solved by DLASDQ.

               dlaset('A', NSIZE, NSIZE, ZERO, ONE, WORK( VT+ST1 ), N );
               dlasdq('U', 0, NSIZE, NSIZE, 0, NRHS, D( ST ), E( ST ), WORK( VT+ST1 ), N, WORK( NWORK ), N, B( ST, 1 ), LDB, WORK( NWORK ), INFO );
               if ( INFO != 0 ) {
                  return;
               }
               dlacpy('A', NSIZE, NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N );
            } else {

               // A large problem. Solve it using divide and conquer.

               dlasda(ICMPQ1, SMLSIZ, NSIZE, SQRE, D( ST ), E( ST ), WORK( U+ST1 ), N, WORK( VT+ST1 ), IWORK( K+ST1 ), WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), WORK( Z+ST1 ), WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ), INFO );
               if ( INFO != 0 ) {
                  return;
               }
               BXST = BX + ST1;
               dlalsa(ICMPQ2, SMLSIZ, NSIZE, NRHS, B( ST, 1 ), LDB, WORK( BXST ), N, WORK( U+ST1 ), N, WORK( VT+ST1 ), IWORK( K+ST1 ), WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), WORK( Z+ST1 ), WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ), INFO );
               if ( INFO != 0 ) {
                  return;
               }
            }
            ST = I + 1;
         }
      } // 60

      // Apply the singular values and treat the tiny ones as zero.

      TOL = RCND*ABS( D( IDAMAX( N, D, 1 ) ) );

      for (I = 1; I <= N; I++) { // 70

         // Some of the elements in D can be negative because 1-by-1
         // subproblems were not solved explicitly.

         if ( ( D( I ) ).abs() <= TOL ) {
            dlaset('A', 1, NRHS, ZERO, ZERO, WORK( BX+I-1 ), N );
         } else {
            RANK = RANK + 1;
            dlascl('G', 0, 0, D( I ), ONE, 1, NRHS, WORK( BX+I-1 ), N, INFO );
         }
         D( I ) = ( D( I ) ).abs();
      } // 70

      // Now apply back the right singular vectors.

      ICMPQ2 = 1;
      for (I = 1; I <= NSUB; I++) { // 80
         ST = IWORK( I );
         ST1 = ST - 1;
         NSIZE = IWORK( SIZEI+I-1 );
         BXST = BX + ST1;
         if ( NSIZE == 1 ) {
            dcopy(NRHS, WORK( BXST ), N, B( ST, 1 ), LDB );
         } else if ( NSIZE <= SMLSIZ ) {
            dgemm('T', 'N', NSIZE, NRHS, NSIZE, ONE, WORK( VT+ST1 ), N, WORK( BXST ), N, ZERO, B( ST, 1 ), LDB );
         } else {
            dlalsa(ICMPQ2, SMLSIZ, NSIZE, NRHS, WORK( BXST ), N, B( ST, 1 ), LDB, WORK( U+ST1 ), N, WORK( VT+ST1 ), IWORK( K+ST1 ), WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), WORK( Z+ST1 ), WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ), INFO );
            if ( INFO != 0 ) {
               return;
            }
         }
      } // 80

      // Unscale and sort the singular values.

      dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO );
      dlasrt('D', N, D, INFO );
      dlascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO );

      return;
      }
