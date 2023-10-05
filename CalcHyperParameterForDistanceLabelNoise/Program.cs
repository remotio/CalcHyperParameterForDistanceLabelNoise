using System;

class program
{



    public static void Main()
    {
        Boundary bd= new Boundary();
        double sigma = 0.25;
        int iteration = 720;
        double goalProb = 0.1;
        decimal expectedValue=bd.Operator(sigma, iteration);
        Console.WriteLine((decimal)goalProb / expectedValue);
    }
}
class Boundary
{
    public int nGrid = 0;
    List<List<double>> gridPoints = new List<List<double>>();
    public int degree = 0;
    List<double> centerCoordinate = new List<double> { 0.5, 0.5 };

    public Boundary() 
    {
        GenerateGridPoint(2);
    }

    public decimal Operator(double sigma,int iteration)
    {
        decimal ans = 0;
        for(int i = 0; i < iteration; i++)
        {
            decimal ave = 0;
            for(int j = 0;j<nGrid;j++)
            {
                List<double>convertedPattern=ConvertCoordinate(gridPoints[j]);
                double dist= Math.Abs(convertedPattern[0] + convertedPattern[1] - 1) / Math.Sqrt(2);
                double probabilityDensityFunction = 1 / (Math.Sqrt(2 * Math.PI) * sigma) * Math.Exp(-dist * dist / (2 * sigma * sigma));
                ave+= (decimal)probabilityDensityFunction;
            }
            ans += ave/nGrid;
            degree++;
        }

        return ans/iteration;
    }


    /// <summary>
    /// Generate list of coordinates of grid points for evaluating model
    /// </summary>
    /// <param name="dim">dimension of the grid points</param>
    private void GenerateGridPoint(int dim)
    {
        int nDivisions = 100;
        nDivisions++;
        nGrid = (int)Math.Pow(nDivisions, dim);

        while (gridPoints.Count < nGrid)
        {
            int count = gridPoints.Count;
            List<double> tmpList = Enumerable.Repeat((double)0, dim).ToList();
            for (int i = 0; i < dim; i++)
            {
                tmpList[i] = count / (int)Math.Pow(nDivisions, dim - i - 1);
                count -= (int)tmpList[i] * (int)Math.Pow(nDivisions, dim - i - 1);
                tmpList[i] /= (double)(nDivisions - 1);
            }
            gridPoints.Add(tmpList);
        }

    }


    /// <summary>
    /// Perform coordinate conversion by rotation matrix (assume 2-dimensional space).
    /// </summary>
    /// <param name="coordinate">coordinate you want to convert</param>
    /// <returns>converted coordinate</returns>
    private List<double> ConvertCoordinate(List<double> coordinate, bool isBoundaryMode = false)
    {
        // Assume 2-dimensional space
        // 境界線ではなく，点を回転させる考え方であるためdegreeに-1をかける
        double radian = -1 * degree * Math.PI / 180.0;
        if (isBoundaryMode)
            radian *= -1;

        List<double> convertedCoord = new List<double>
            {
                (coordinate[0] - centerCoordinate[0]) * Math.Cos(radian) - (coordinate[1] - centerCoordinate[1]) * Math.Sin(radian) + centerCoordinate[0],
                (coordinate[0] - centerCoordinate[0]) * Math.Sin(radian) + (coordinate[1] - centerCoordinate[1]) * Math.Cos(radian) + centerCoordinate[1]
            };

        return convertedCoord;
    }
}

