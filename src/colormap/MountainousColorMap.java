package colormap;
import java.awt.Color;


public class MountainousColorMap implements Mapper<Float, Color>
{
	@Override
	public Color map(Float value)
	{
		//Ocean   .7f .8f  .3f
		//Beach   .2f .85f .85f
		//Grass: .28f  .8f .7f
		//Brown: .09, .85f .8f
		//Snow:   DC,  0f,  1f
		
		float hue, saturation, brightness;
		
		if (value < 1f/2f)
		{
			value = value * 2;
			value = value * value * (3f - 2f * value);
			value = (float)Math.pow(value, 3);
			
			hue        =  .2f * value + (1f - value) * .65f;
			saturation =  .4f * value + (1f - value) *  .8f;
			brightness =  .9f * value + (1f - value) *  .3f;
		}
		else if (value < 4f/6f)
		{
			value = (float)Math.pow(value * 6f - 3f, .7f);
			
			hue        = .28f * value + (1f - value) *  .2f;
			saturation =  .8f * value + (1f - value) * .85f;
			brightness =  .7f * value + (1f - value) * .85f;
		}
		else if (value < 5f/6f)
		{
			value = (float)Math.pow(value * 6f - 4f, 2.5f);
			
			hue        = .09f * value + (1f - value) * .28f;
			saturation = .85f * value + (1f - value) *  .8f;
			brightness =  .7f * value + (1f - value) *  .7f;
		}
		else
		{
			value = (float)Math.pow(value * 6f - 5f, 2f);
			
			hue        = .07f * value + (1f - value) * .09f;
			saturation =   0f * value + (1f - value) *  .85f;
			brightness =   1f * value + (1f - value) *  .7f;
		}
		
		return new Color(Color.HSBtoRGB(hue, saturation, brightness));
	}

}
