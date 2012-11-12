package squarediamond;

import java.util.Random;

import src.Array2D;


public final class AveragedSquareDiamondStyle implements SquareDiamondArray2DPopulator.ISquareDiamondStyle
{
	private final Random random;
	
	public AveragedSquareDiamondStyle(Random random)
	{
		this.random = random;
	}

	@Override
	public void setup(Array2D values)
	{
	}
	
	@Override
	public void onNextIteration()
	{
	}
	
	@Override
	public float getSeedValue(int i, int j)
	{
		return random.nextFloat();
	}

	@Override
	public float squareValue(float topLeft, float topRight, float bottomLeft, float bottomRight)
	{
		return (topLeft + topRight + bottomLeft + bottomRight) * 1f/4f;
	}

	@Override
	public float diamondValue(float left, float right, float top, float bottom)
	{
		return (left + right + top + bottom) * 1f/4f;
	}

	@Override
	public float diamondValue_top(float left, float right, float bottom)
	{
		return (left + right + bottom) * 1f/3f;
	}

	@Override
	public float diamondValue_bottom(float left, float right, float top)
	{
		return (left + right + top) * 1f/3f;
	}

	@Override
	public float diamondValue_left(float right, float top, float bottom)
	{
		return (right + top + bottom) * 1f/3f;
	}

	@Override
	public float diamondValue_right(float left, float top, float bottom)
	{
		return (left + top + bottom) * 1f/3f;
	}

}
