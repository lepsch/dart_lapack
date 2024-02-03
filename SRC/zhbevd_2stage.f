      SUBROUTINE ZHBEVD_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             RWORK( * ), W( * );
      COMPLEX*16         AB( LDAB, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D0, 0.0D0 ), CONE = ( 1.0D0, 0.0D0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, IMAX, INDE, INDWK2, INDRWK, ISCALE, LLWORK, INDWK, LHTRD, LWTRD, IB, INDHOUS, LIWMIN, LLRWK, LLWK2, LRWMIN, LWMIN;
      double             ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      double             DLAMCH, ZLANHB;
      // EXTERNAL LSAME, DLAMCH, ZLANHB, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSTERF, XERBLA, ZGEMM, ZLACPY, ZLASCL, ZSTEDC, ZHETRD_HB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 .OR. LRWORK.EQ.-1 )

      INFO = 0
      if ( N.LE.1 ) {
         LWMIN = 1
         LRWMIN = 1
         LIWMIN = 1
      } else {
         IB    = ILAENV2STAGE( 2, 'ZHETRD_HB2ST', JOBZ, N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'ZHETRD_HB2ST', JOBZ, N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'ZHETRD_HB2ST', JOBZ, N, KD, IB, -1 )
         if ( WANTZ ) {
            LWMIN = 2*N**2
            LRWMIN = 1 + 5*N + 2*N**2
            LIWMIN = 3 + 5*N
         } else {
            LWMIN  = MAX( N, LHTRD + LWTRD )
            LRWMIN = N
            LIWMIN = 1
         }
      }
      if ( .NOT.( LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( KD.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -6
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -9
      }

      if ( INFO.EQ.0 ) {
         WORK( 1 )  = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -11
         } else if ( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) {
            INFO = -13
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -15
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('ZHBEVD_2STAGE', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;

      if ( N.EQ.1 ) {
         W( 1 ) = DBLE( AB( 1, 1 ) )
         if (WANTZ) Z( 1, 1 ) = CONE;
         RETURN
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS    = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = ZLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
      ISCALE = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE.EQ.1 ) {
         if ( LOWER ) {
            zlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            zlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
      }

      // Call ZHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.

      INDE    = 1
      INDRWK  = INDE + N
      LLRWK   = LRWORK - INDRWK + 1
      INDHOUS = 1
      INDWK   = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWK + 1
      INDWK2  = INDWK + N*N
      LLWK2   = LWORK - INDWK2 + 1

      zhetrd_hb2st("N", JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWK ), LLWORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEDC.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, RWORK( INDE ), INFO );
      } else {
         zstedc('I', N, W, RWORK( INDE ), WORK, N, WORK( INDWK2 ), LLWK2, RWORK( INDRWK ), LLRWK, IWORK, LIWORK, INFO );
         zgemm('N', 'N', N, N, N, CONE, Z, LDZ, WORK, N, CZERO, WORK( INDWK2 ), N );
         zlacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = N
         } else {
            IMAX = INFO - 1
         }
         dscal(IMAX, ONE / SIGMA, W, 1 );
      }

      WORK( 1 )  = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of ZHBEVD_2STAGE

      }
