package squarediamond;
import java.util.Random;

import squarediamond.SquareDiamondArray2DPopulator.ISquareDiamondStyle;
import src.Array2D;


public abstract class InterpolatedSquareDiamondStyle implements ISquareDiamondStyle
{
	protected final Random random;
	
	public InterpolatedSquareDiamondStyle(Random random)
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
	public final float getSeedValue(int i, int j)
	{
		return random.nextFloat();
	}
	
	protected abstract float nextWeight();
	
	@Override
	public final float squareValue(float topLeft, float topRight, float bottomLeft, float bottomRight)
	{
		float weight;
		
		float left  = topLeft  * (weight = nextWeight()) + (1f - weight) * bottomLeft;
		float right = topRight * (weight = nextWeight()) + (1f - weight) * bottomRight;
		
		return left * (weight = nextWeight()) + (1f - weight) * right;
	}
	
	@Override
	public final float diamondValue(float left, float right, float top, float bottom)
	{
		float weight;
		float a = left * (weight = nextWeight()) + (1f - weight) * right;
		float b = top * (weight = nextWeight()) + (1f - weight) * bottom;
		
		return a * (weight = nextWeight()) + (1f - weight) * b;
	}

	@Override
	public final float diamondValue_top(float left, float right, float bottom)
	{
		float weight;
		float a = left * (weight = nextWeight()) + (1f - weight) * right;
		
		return a * (weight = nextWeight()) + (1f - weight) * bottom;
	}

	@Override
	public final float diamondValue_bottom(float left, float right, float top)
	{
		float weight;
		float a = left * (weight = nextWeight()) + (1f - weight) * right;
		
		return a * (weight = nextWeight()) + (1f - weight) * top;
	}

	@Override
	public final float diamondValue_left(float right, float top, float bottom)
	{
		float weight;
		float b = top * (weight = nextWeight()) + (1f - weight) * bottom;
		
		return right * (weight = nextWeight()) + (1f - weight) * b;
	}

	@Override
	public final float diamondValue_right(float left, float top, float bottom)
	{
		float weight;
		float b = top * (weight = nextWeight()) + (1f - weight) * bottom;
		
		return left * (weight = nextWeight()) + (1f - weight) * b;
	}

}
