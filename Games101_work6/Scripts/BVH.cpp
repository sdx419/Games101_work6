#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    root = recursiveBuildWithSAH(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim) {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x <
                       f2->getBounds().Centroid().x;
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y <
                       f2->getBounds().Centroid().y;
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z <
                       f2->getBounds().Centroid().z;
            });
            break;
        }

        auto beginning = objects.begin();
        auto middling = objects.begin() + (objects.size() / 2);
        auto ending = objects.end();

        auto leftshapes = std::vector<Object*>(beginning, middling);
        auto rightshapes = std::vector<Object*>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }

    return node;
}

BVHBuildNode* BVHAccel::recursiveBuildWithSAH(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode;

    if (objects.size() == 1)
    {
        node->bounds = objects[0]->getBounds();
        node->left = nullptr;
        node->right = nullptr;
        node->object = objects[0];
        return node;
    }
    else if (objects.size() == 2)
    {
        node->left = recursiveBuildWithSAH(std::vector{ objects[0] });
        node->right = recursiveBuildWithSAH(std::vector{ objects[1] });
        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else
    {
        Bounds3 bounds;
        for (int i = 0; i < objects.size(); i++)
            bounds = Union(bounds, objects[i]->getBounds());
        node->bounds = bounds;

        int dim = bounds.maxExtent();
        switch (dim)
        {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                Bounds3 b1 = f1->getBounds();
                Bounds3 b2 = f2->getBounds();
                return b1.pMax.x < b2.pMax.x;
                });
            break;

        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                Bounds3 b1 = f1->getBounds();
                Bounds3 b2 = f2->getBounds();
                return b1.pMax.y < b2.pMax.y;
                });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                Bounds3 b1 = f1->getBounds();
                Bounds3 b2 = f2->getBounds();
                return b1.pMax.z < b2.pMax.z;
                });
            break;
            
        default:
            break;
        }

        int partIndex = 0;
        auto surfaceArea = std::numeric_limits<double>::max();

        for (int partitionIndex = 0; partitionIndex < objects.size() - 1; partitionIndex++)
        {
            Bounds3 bound1, bound2;

            for (int i = 0; i <= partitionIndex; i++)
                bound1 = Union(bound1, objects[i]->getBounds());

            for (int i = partitionIndex + 1; i < objects.size(); i++)
                bound2 = Union(bound2, objects[i]->getBounds());

            if (bound1.SurfaceArea() / bounds.SurfaceArea() * partitionIndex + bound2.SurfaceArea() / bounds.SurfaceArea() * (objects.size() - partitionIndex) < surfaceArea)
            {
                partIndex = partitionIndex;
                surfaceArea = bound1.SurfaceArea() / bounds.SurfaceArea() * partitionIndex + bound2.SurfaceArea() / bounds.SurfaceArea() * (objects.size() - partitionIndex);
            }
        }

        auto leftShapes = std::vector<Object*>(objects.begin(), objects.begin() + partIndex + 1);
        auto rightShapes = std::vector<Object*>(objects.begin() + partIndex + 1, objects.end());
        node->left = recursiveBuildWithSAH(leftShapes);
        node->right = recursiveBuildWithSAH(rightShapes);
        return node;
    }
    



}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const  
{  
   // TODO Traverse the BVH to find intersection  
    std::array<int, 3> dirIsNeg = { 1,1,1 };
    if (ray.direction.x < 0)
        dirIsNeg[0] = 0;

    if (ray.origin.y < 0)
        dirIsNeg[1] = 0;

    if (ray.origin.z > 0)
        dirIsNeg[2] = 0;

   if (!node->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg))  
   {  
       Intersection inter;  
       inter.happened = false;  
       return inter;  
   }  

   if (node->left == nullptr && node->right == nullptr && node->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg))  
   {  
       //std::cout << "intersection occured\n";
       return node->object->getIntersection(ray);  
   }  

   Intersection leftIntersection, rightIntersection;  

   if (node->left != nullptr && node->left->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg))  
   {  
       leftIntersection = getIntersection(node->left, ray);
   }  

   if (node->right != nullptr && node->right->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg))
   {
       rightIntersection = getIntersection(node->right, ray);
   }

   return (leftIntersection.happened && (!rightIntersection.happened || leftIntersection.distance < rightIntersection.distance)) 
       ? leftIntersection : rightIntersection;
}