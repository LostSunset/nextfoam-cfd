// Helper: quantize point to integer key for hashing
struct PointHash
{
    scalar tol;
    PointHash(scalar t) : tol(max(t, SMALL)) {}

    std::size_t operator()(const point& p) const noexcept
    {
        scalar invTol = 1.0 / tol;
        // avoid NaN propagation
        scalar px = std::isfinite(p.x()) ? p.x() : 0.0;
        scalar py = std::isfinite(p.y()) ? p.y() : 0.0;
        scalar pz = std::isfinite(p.z()) ? p.z() : 0.0;

        long long ix = llround(px * invTol);
        long long iy = llround(py * invTol);
        long long iz = llround(pz * invTol);

        // standard hash combination
        std::size_t h = 1469598103934665603ULL;
        h ^= ix + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= iy + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= iz + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};


struct PointEq
{
    scalar tol;
    PointEq(scalar t) : tol(max(t, SMALL)) {}
    bool operator()(const point& a, const point& b) const noexcept
    {
        return mag(a - b) < tol;
    }
};


bool coincidentGeometry
(
    const polyMesh& mesh,
    const label patchA,
    const label patchB,
    const scalar tol = 1e-6
)
{
    const polyPatch& pA = mesh.boundaryMesh()[patchA];
    const polyPatch& pB = mesh.boundaryMesh()[patchB];

    const pointField& ptsA = pA.points();
    const pointField& ptsB = pB.points();

    // 1. Bounding box check
    boundBox boxA(ptsA);
    boundBox boxB(ptsB);

    if (!boxA.overlaps(boxB))
    {
        return false;
    }

    // 2. Area + centroid
    scalar areaA = 0, areaB = 0;
    vector centroidA(Zero), centroidB(Zero);

    forAll(pA, i)
    {
        scalar a = pA[i].mag(mesh.points());
        areaA += a;
        centroidA += pA[i].centre(mesh.points()) * a;
    }
    centroidA /= (areaA + SMALL);

    forAll(pB, i)
    {
        scalar a = pB[i].mag(mesh.points());
        areaB += a;
        centroidB += pB[i].centre(mesh.points()) * a;
    }
    centroidB /= (areaB + SMALL);

    if (mag(areaA - areaB) > tol*areaA)
    {
        return false;
    }

    if (mag(centroidA - centroidB) > tol)
    {
        return false;
    }

/*    // 3. Build hash set for patch B points
    std::unordered_set<point, PointHash, PointEq> setB
    (
        0,
        PointHash(tol),
        PointEq(tol)
    );

    forAll(ptsB, j)
    {
        const point& p = ptsB[j];
        if (finite(p.x()) && finite(p.y()) && finite(p.z()))
        {
            setB.insert(p);
        }
    }

    // 4. Check each A point has match in B
    forAll(ptsA, i)
    {
        const point& p = ptsA[i];
        if (finite(p.x()) && finite(p.y()) && finite(p.z()))
        {
            if (setB.find(p) == setB.end())
                return false;
        }
    } */

    return true;
}

