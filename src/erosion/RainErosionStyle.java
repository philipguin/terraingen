package erosion;

import java.util.Random;

import src.Array2D;
import src.IDerivator;

public class RainErosionStyle implements IErosionStyle
{
	private final Point point = new Point(0, 0);
	
	private final Random random;
	private final IDerivator<Array2D> distributionMaker;
	private final float initialVelocity;
	
	private Array2D distribution;
	
	public RainErosionStyle(Random random, IDerivator<Array2D> distributionMaker, float initialVelocity)
	{
		this.random = random;
		this.distributionMaker = distributionMaker;
		this.initialVelocity = initialVelocity;
	}
	
	@Override
	public void setup(Array2D heightMap)
	{
		distribution = distributionMaker.derive(heightMap);
	}

	@Override
	public Point makeStartPoint()
	{
		do
		{
			point.x = random.nextInt(distribution.width - 1);
			point.y = random.nextInt(distribution.height - 1);
		}
		while (random.nextFloat() >= distribution.get(point.x, point.y));
		
		return point;
	}

	@Override
	public float makeInitialVelocity()
	{
		return initialVelocity;
	}

	@Override
	public float makeInitialSediment()
	{
		return 0;
	}

}
