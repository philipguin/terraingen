package erosion;

import java.util.Random;

import src.Array2D;

public class DropletErosionStyle implements IErosionStyle
{
	private final Point point = new Point(0, 0);
	
	private final Random random;
	private int width, height;
	private final float initialVelocity, initialSediment;
	
	public DropletErosionStyle(Random random, float initialVelocity, float initialSediment)
	{
		this.random = random;
		this.initialVelocity = initialVelocity;
		this.initialSediment = initialSediment;
	}
	
	@Override
	public void setup(Array2D heightMap)
	{
		width = heightMap.width;
		height = heightMap.height;
	}

	@Override
	public Point makeStartPoint()
	{
		point.x = random.nextInt(width - 1);
		point.y = random.nextInt(height - 1);
		
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
		return initialSediment;
	}

}
