      SUBROUTINE SSBEVD_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LQUERY, WANTZ;
      int                IINFO, INDE, INDWK2, INDWRK, ISCALE, LIWMIN, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS, LLWRK2;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SLAMCH, SLANSB, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, SLANSB, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLASCL, SSCAL, SSTEDC, SSTERF, XERBLA, SSYTRD_SB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

      INFO = 0
      if ( N.LE.1 ) {
         LIWMIN = 1
         LWMIN = 1
      } else {
         IB    = ILAENV2STAGE( 2, 'SSYTRD_SB2ST', JOBZ, N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'SSYTRD_SB2ST', JOBZ, N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'SSYTRD_SB2ST', JOBZ, N, KD, IB, -1 )
         if ( WANTZ ) {
            LIWMIN = 3 + 5*N
            LWMIN = 1 + 5*N + 2*N**2
         } else {
            LIWMIN = 1
            LWMIN = MAX( 2*N, N+LHTRD+LWTRD )
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
         WORK( 1 )  = SROUNDUP_LWORK(LWMIN)
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -11
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -13
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('SSBEVD_2STAGE', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;

      if ( N.EQ.1 ) {
         W( 1 ) = AB( 1, 1 )
         if (WANTZ) Z( 1, 1 ) = ONE;
         RETURN
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS    = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = SLANSB( 'M', UPLO, N, KD, AB, LDAB, WORK )
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
            slascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            slascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
      }

      // Call SSYTRD_SB2ST to reduce band symmetric matrix to tridiagonal form.

      INDE    = 1
      INDHOUS = INDE + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1
      INDWK2  = INDWRK + N*N
      LLWRK2  = LWORK - INDWK2 + 1

      ssytrd_sb2st("N", JOBZ, UPLO, N, KD, AB, LDAB, W, WORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, call SSTEDC.

      if ( .NOT.WANTZ ) {
         ssterf(N, W, WORK( INDE ), INFO );
      } else {
         sstedc('I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO )          CALL SGEMM( 'N', 'N', N, N, N, ONE, Z, LDZ, WORK( INDWRK ), N, ZERO, WORK( INDWK2 ), N );
         slacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      if (ISCALE.EQ.1) CALL SSCAL( N, ONE / SIGMA, W, 1 );

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of SSBEVD_2STAGE

      }
