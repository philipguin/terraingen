package squarediamond;

import squarediamond.SquareDiamondArray2DPopulator.ISquareDiamondStyle;
import src.Array2D;


public abstract class SimpleCompositeSquareDiamondStyle implements ISquareDiamondStyle
{
	private final ISquareDiamondStyle base;
	
	public SimpleCompositeSquareDiamondStyle(ISquareDiamondStyle base)
	{
		this.base = base;
	}
	
	@Override
	public void setup(Array2D values)
	{
		base.setup(values);
	}
	
	@Override
	public void onNextIteration()
	{
		base.onNextIteration();
	}

	@Override
	public final float getSeedValue(int i, int j)
	{
		return base.getSeedValue(i, j);
	}

	protected abstract float bias(float value);
	
	@Override
	public final float squareValue(float topLeft, float topRight, float bottomLeft, float bottomRight)
	{
		return bias(base.squareValue(topLeft, topRight, bottomLeft, bottomRight));
	}

	@Override
	public final float diamondValue(float left, float right, float top, float bottom)
	{
		return bias(base.diamondValue(left, right, top, bottom));
	}

	@Override
	public final float diamondValue_top(float left, float right, float bottom)
	{
		return bias(base.diamondValue_top(left, right, bottom));
	}

	@Override
	public final float diamondValue_bottom(float left, float right, float top)
	{
		return bias(base.diamondValue_bottom(left, right, top));
	}

	@Override
	public final float diamondValue_left(float right, float top, float bottom)
	{
		return bias(base.diamondValue_left(right, top, bottom));
	}

	@Override
	public final float diamondValue_right(float left, float top, float bottom)
	{
		return bias(base.diamondValue_right(left, top, bottom));
	}

}
