      import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlaqz0.dart';

void zlaqz2( final bool ILSCHUR, final bool ILQ, final bool ILZ, final int N, final int ILO, final int IHI, final int NW,
      final Matrix<Complex> A_, final int LDA, final Matrix<Complex> B_, final int LDB, final Matrix<Complex> Q_, final int LDQ, final Matrix<Complex> Z_, final int LDZ,
      final Box<int> NS, final Box<int> ND, final Array<Complex> ALPHA_, final Array<Complex> BETA_, final Matrix<Complex> QC_, final int LDQC, final Matrix<Complex> ZC_, final int LDZC,
      final Array<Complex> WORK_, final int LWORK, final Array<double> RWORK_, final int REC, final Box<int> INFO, ){
final A=A_.dim(LDA);
final B=B_.dim(LDB);
final Q=Q_.dim(LDQ);
final Z=Z_.dim(LDZ);
final ALPHA=ALPHA_.dim();
final BETA=BETA_.dim();
final QC=QC_.dim(LDQC);
final ZC=ZC_.dim(LDZC);
final WORK=WORK_.dim();
final RWORK=RWORK_.dim();
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local Scalars
      int      JW, KWTOP, KWBOT, ISTOPM, ISTARTM, K, K2, ZTGEXC_INFO, IFST, ILST, LWORKREQ, QZ_SMALL_INFO;
      double           SMLNUM, ULP, SAFMIN, SAFMAX, C1, TEMPR;
      Complex  S, S1, TEMP;

      // External Functions
      // EXTERNAL :: XERBLA, ZLAQZ0, ZLAQZ1, ZLACPY, ZLASET, ZGEMM, ZTGEXC, ZLARTG, ZROT
      // double          , EXTERNAL :: DLAMCH;

      INFO.value = 0;

      // Set up deflation window
      JW = min( NW, IHI-ILO+1 );
      KWTOP = IHI-JW+1;
      if ( KWTOP == ILO ) {
         S = Complex.zero;
      } else {
         S = A[ KWTOP][KWTOP-1 ];
      }

      // Determine required workspace
      IFST = 1;
      ILST = JW;
      ZLAQZ0('S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHA, BETA, QC, LDQC, ZC, LDZC, WORK, -1, RWORK, REC+1, QZ_SMALL_INFO );
      LWORKREQ = INT( WORK( 1 ) )+2*JW**2;
      LWORKREQ = max( LWORKREQ, N*NW, 2*NW**2+N );
      if ( LWORK == -1 ) {
         // workspace query, quick return;
         WORK[1] = LWORKREQ;
         return;
      } else if ( LWORK < LWORKREQ ) {
         INFO.value = -26;
      }

      if ( INFO.value != 0 ) {
         xerbla('ZLAQZ2', -INFO.value );
         return;
      }

      // Get machine constants
      SAFMIN = dlamch( 'SAFE MINIMUM' );
      SAFMAX = ONE/SAFMIN;
      ULP = dlamch( 'PRECISION' );
      SMLNUM = SAFMIN*( N.toDouble()/ULP );

      if ( IHI == KWTOP ) {
         // 1 by 1 deflation window, just try a regular deflation
         ALPHA[KWTOP] = A( KWTOP, KWTOP );
         BETA[KWTOP] = B( KWTOP, KWTOP );
         NS.value = 1;
         ND.value = 0;
         if ( ( S ).abs() <= max( SMLNUM, ULP*( A( KWTOP, KWTOP ) ).abs() ) ) {
            NS.value = 0;
            ND.value = 1;
            if ( KWTOP > ILO ) {
               A[KWTOP][KWTOP-1] = Complex.zero;
            }
         }
      }


      // Store window in case of convergence failure
      zlacpy('ALL', JW, JW, A( KWTOP, KWTOP ), LDA, WORK, JW );
      zlacpy('ALL', JW, JW, B( KWTOP, KWTOP ), LDB, WORK( JW**2+ 1 ), JW );

      // Transform window to real schur form
      zlaset('FULL', JW, JW, Complex.zero, Complex.one, QC, LDQC );
      zlaset('FULL', JW, JW, Complex.zero, Complex.one, ZC, LDZC );
      zlaqz0('S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHA, BETA, QC, LDQC, ZC, LDZC, WORK( 2*JW**2+1 ), LWORK-2*JW**2, RWORK, REC+1, QZ_SMALL_INFO );

      if ( QZ_SMALL_INFO != 0 ) {
         // Convergence failure, restore the window and exit
         ND.value = 0;
         NS.value = JW-QZ_SMALL_INFO;
         zlacpy('ALL', JW, JW, WORK, JW, A( KWTOP, KWTOP ), LDA );
         zlacpy('ALL', JW, JW, WORK( JW**2+1 ), JW, B( KWTOP, KWTOP ), LDB );
         return;
      }

      // Deflation detection loop
      if ( KWTOP == ILO || S == Complex.zero ) {
         KWBOT = KWTOP-1;
      } else {
         KWBOT = IHI;
         K = 1;
         K2 = 1;
         while (K <= JW) {
               // Try to deflate eigenvalue
               TEMPR = ( A( KWBOT, KWBOT ) ).abs();
               if ( TEMPR == ZERO ) {
                  TEMPR = ( S ).abs();
               }
               if ( ( ( S*QC( 1, KWBOT-KWTOP+1 ) ).abs() ) <= max( ULP* TEMPR, SMLNUM ) ) {
                  // Deflatable
                  KWBOT = KWBOT-1;
               } else {
                  // Not deflatable, move out of the way
                  IFST = KWBOT-KWTOP+1;
                  ILST = K2;
                  ztgexc( true , true , JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, QC, LDQC, ZC, LDZC, IFST, ILST, ZTGEXC_INFO );
                  K2 = K2+1;
               }

               K = K+1;
         }
      }

      // Store eigenvalues
      ND.value = IHI-KWBOT;
      NS.value = JW-ND.value;
      K = KWTOP;
      while (K <= IHI) {
         ALPHA[K] = A( K, K );
         BETA[K] = B( K, K );
         K = K+1;
      }

      if ( KWTOP != ILO && S != Complex.zero ) {
         // Reflect spike back, this will create optimally packed bulges
         A[KWTOP:KWBOT][KWTOP-1] = A( KWTOP, KWTOP-1 ) *DCONJG( QC( 1, 1:JW-ND.value ) );
         for (K = KWBOT-1; K >= KWTOP; K--) {
            zlartg(A( K, KWTOP-1 ), A( K+1, KWTOP-1 ), C1, S1, TEMP );
            A[K][KWTOP-1] = TEMP;
            A[K+1][KWTOP-1] = Complex.zero;
            K2 = max( KWTOP, K-1 );
            zrot(IHI-K2+1, A( K, K2 ), LDA, A( K+1, K2 ), LDA, C1, S1 );
            zrot(IHI-( K-1 )+1, B( K, K-1 ), LDB, B( K+1, K-1 ), LDB, C1, S1 );
            zrot(JW, QC( 1, K-KWTOP+1 ), 1, QC( 1, K+1-KWTOP+1 ), 1, C1, DCONJG( S1 ) );
         }

         // Chase bulges down
         ISTARTM = KWTOP;
         ISTOPM = IHI;
         K = KWBOT-1;
         while (K >= KWTOP) {

            // Move bulge down and remove it
            for (K2 = K; K2 <= KWBOT-1; K2++) {
               zlaqz1( true , true , K2, KWTOP, KWTOP+JW-1, KWBOT, A, LDA, B, LDB, JW, KWTOP, QC, LDQC, JW, KWTOP, ZC, LDZC );
            }

            K = K-1;
         }

      }

      // Apply Qc and Zc to rest of the matrix
      if ( ILSCHUR ) {
         ISTARTM = 1;
         ISTOPM = N;
      } else {
         ISTARTM = ILO;
         ISTOPM = IHI;
      }

      if ( ISTOPM-IHI > 0 ) {
         zgemm('C', 'N', JW, ISTOPM-IHI, JW, Complex.one, QC, LDQC, A( KWTOP, IHI+1 ), LDA, Complex.zero, WORK, JW );
         zlacpy('ALL', JW, ISTOPM-IHI, WORK, JW, A( KWTOP, IHI+1 ), LDA );
         zgemm('C', 'N', JW, ISTOPM-IHI, JW, Complex.one, QC, LDQC, B( KWTOP, IHI+1 ), LDB, Complex.zero, WORK, JW );
         zlacpy('ALL', JW, ISTOPM-IHI, WORK, JW, B( KWTOP, IHI+1 ), LDB );
      }
      if ( ILQ ) {
         zgemm('N', 'N', N, JW, JW, Complex.one, Q( 1, KWTOP ), LDQ, QC, LDQC, Complex.zero, WORK, N );
         zlacpy('ALL', N, JW, WORK, N, Q( 1, KWTOP ), LDQ );
      }

      if ( KWTOP-1-ISTARTM+1 > 0 ) {
         zgemm('N', 'N', KWTOP-ISTARTM, JW, JW, Complex.one, A( ISTARTM, KWTOP ), LDA, ZC, LDZC, Complex.zero, WORK, KWTOP-ISTARTM );
        zlacpy('ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, A( ISTARTM, KWTOP ), LDA );
         zgemm('N', 'N', KWTOP-ISTARTM, JW, JW, Complex.one, B( ISTARTM, KWTOP ), LDB, ZC, LDZC, Complex.zero, WORK, KWTOP-ISTARTM );
        zlacpy('ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, B( ISTARTM, KWTOP ), LDB );
      }
      if ( ILZ ) {
         zgemm('N', 'N', N, JW, JW, Complex.one, Z( 1, KWTOP ), LDZ, ZC, LDZC, Complex.zero, WORK, N );
         zlacpy('ALL', N, JW, WORK, N, Z( 1, KWTOP ), LDZ );
      }

}
